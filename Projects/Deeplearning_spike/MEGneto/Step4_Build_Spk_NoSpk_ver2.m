%% ------------------------------------------------------------------------
% Step4_Build_Spike_NoSpike_dataset.m
%
% Build Spike/NoSpike DL-ready dataset from per-reviewer event MAT files
% + meta txt files.
%
% Input folder structure:
%   root/Adi/*_spike_events.mat
%   root/Adi/*_nospike_events.mat
%   root/Adi/meta_files/*_Adi_meta_data.txt
%   ...
%
% Output:
%   root/dataset_DL_ready.mat
%   root/Step4_file_read_summary.csv
%   root/Step4_dataset_summary.csv
%
% Event MAT format expected:
%   events: N x time x channels
%
% Final dataset format:
%   dataset.X_train: N x channels x time
% -------------------------------------------------------------------------

clc; clear;

%% ======================== USER SETTINGS =================================

opts.root = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files';

opts.raters = {'Adi','Josh','Manoj','Pradeep'};

% opts.outFile = fullfile(opts.root, 'dataset_DL_ready.mat');
opts.outFile = fullfile(opts.root, 'dataset_DL_ready_win_m250_p500.mat');

% Expected epoch dimensions
opts.cropSamp = 1:306;        % time samples to keep
% opts.expectedNTime = 306;
%%
Fs = 200;
origWinSec = [-0.500 1.025];
desiredWinSec = [-0.250 0.500];

tOrig = origWinSec(1):1/Fs:origWinSec(2);

opts.cropSamp = find(tOrig >= desiredWinSec(1) & tOrig <= desiredWinSec(2));
opts.expectedNTime = numel(opts.cropSamp);
opts.expectedNChan = 306;

fprintf('Using crop window %.3f to %.3f sec: %d samples\n', ...
    desiredWinSec(1), desiredWinSec(2), opts.expectedNTime);

%%
opts.expectedNChan = 306;

opts.Fs_default = 200;
opts.useSingle = true;

% Consensus
opts.minRatersForConsensus = 2;
opts.consensusRule = 'majority';   % ties are dropped

% No-spike handling
opts.addNoSpikeEpochs = true;

% Important for reconstructed data:
% If the same no-spike epoch was saved under two reviewers, do not duplicate it.
opts.deduplicateNoSpike = true;

% Split
opts.rngSeed = 1;
opts.splitRatio = [0.70 0.15 0.15];   % train / val / test

% Plot one example
opts.makeQuickPlot = true;

rng(opts.rngSeed);

%% ======================== CHECK INPUT FOLDERS ============================

if ~exist(opts.root, 'dir')
    error('Root folder not found: %s', opts.root);
end

for r = 1:numel(opts.raters)
    raterDir = fullfile(opts.root, opts.raters{r});
    metaDir  = fullfile(raterDir, 'meta_files');

    if ~exist(raterDir, 'dir')
        warning('Missing rater folder: %s', raterDir);
    end

    if ~exist(metaDir, 'dir')
        warning('Missing meta_files folder: %s', metaDir);
    end
end

%% ======================== INITIALIZE STORAGE =============================

% Positive spike candidates pooled by base + sample + code.
key2idx = containers.Map('KeyType','char','ValueType','double');

X_cell    = {};   % each cell is C x T
subj_cell = {};
base_cell = {};
sample_v  = [];
code_v    = [];
Fs_v      = [];

labelsMat = NaN(0, numel(opts.raters));  % Ncandidate x Nrater

% No-spike extras
X_ns_cell = {};
subj_ns   = {};
base_ns   = {};
Fs_ns     = [];

% For no-spike deduplication
nospikeSeen = containers.Map('KeyType','char','ValueType','logical');

Rows = {};

%% ======================== LOAD PER-RATER FILES ===========================

for r = 1:numel(opts.raters)

    rname = opts.raters{r};
    raterDir = fullfile(opts.root, rname);
    metaDir  = fullfile(raterDir, 'meta_files');

    spkFiles = dir(fullfile(raterDir, '*_spike_events.mat'));

    fprintf('\n==============================\n');
    fprintf('Rater: %s\n', rname);
    fprintf('Spike files found: %d\n', numel(spkFiles));
    fprintf('==============================\n');

    for k = 1:numel(spkFiles)

        spkFn = fullfile(spkFiles(k).folder, spkFiles(k).name);
        base  = erase(spkFiles(k).name, '_spike_events.mat');

        nosFn  = fullfile(raterDir, [base '_nospike_events.mat']);
        metaFn = fullfile(metaDir, sprintf('%s_%s_meta_data.txt', base, rname));

        status = 'OK';
        nSpikeUsed = 0;
        nNoUsed = 0;
        nAgree = NaN;
        nMeta = NaN;

        fprintf('\n--- %s | %s ---\n', rname, base);

        %% ---------- Check required files ----------
        if ~isfile(metaFn)
            warning('Missing meta file: %s', metaFn);
            status = 'MISSING_META';
            Rows(end+1,:) = {rname, base, spkFn, metaFn, nosFn, status, ...
                nSpikeUsed, nNoUsed, nMeta, nAgree}; %#ok<SAGROW>
            continue;
        end

        %% ---------- Load spike MAT ----------
        try
            [Es, metas] = load_event_mat_flexible(spkFn);
        catch ME
            warning('Could not load spike MAT: %s\n%s', spkFn, ME.message);
            status = 'BAD_SPIKE_MAT';
            Rows(end+1,:) = {rname, base, spkFn, metaFn, nosFn, status, ...
                nSpikeUsed, nNoUsed, nMeta, nAgree}; %#ok<SAGROW>
            continue;
        end

        Fs_s = getFs(metas, opts.Fs_default);

        if ndims(Es) ~= 3
            warning('Spike events are not 3D: %s', spkFn);
            status = 'BAD_SPIKE_DIM';
            Rows(end+1,:) = {rname, base, spkFn, metaFn, nosFn, status, ...
                nSpikeUsed, nNoUsed, nMeta, nAgree}; %#ok<SAGROW>
            continue;
        end

        % Es is expected as N x time x channels
        cropIdx = opts.cropSamp(opts.cropSamp <= size(Es,2));
        Es = Es(:, cropIdx, :);

        % Convert to N x channels x time
        Xs = permute(Es, [1 3 2]);

        if size(Xs,2) ~= opts.expectedNChan || size(Xs,3) ~= opts.expectedNTime
            warning('Unexpected spike size for %s: Xs = [%d x %d x %d]', ...
                base, size(Xs,1), size(Xs,2), size(Xs,3));
            status = 'BAD_SPIKE_SIZE';
            Rows(end+1,:) = {rname, base, spkFn, metaFn, nosFn, status, ...
                nSpikeUsed, nNoUsed, nMeta, nAgree}; %#ok<SAGROW>
            continue;
        end

        %% ---------- Read meta ----------
        try
            Tm = read_meta_flexible(metaFn);
        catch ME
            warning('Could not read meta file: %s\n%s', metaFn, ME.message);
            status = 'BAD_META';
            Rows(end+1,:) = {rname, base, spkFn, metaFn, nosFn, status, ...
                nSpikeUsed, nNoUsed, nMeta, nAgree}; %#ok<SAGROW>
            continue;
        end

        requiredVars = {'sample','agree','code'};
        if isempty(Tm) || ~all(ismember(requiredVars, Tm.Properties.VariableNames))
            warning('Meta missing sample/agree/code: %s', metaFn);
            status = 'BAD_META_COLUMNS';
            Rows(end+1,:) = {rname, base, spkFn, metaFn, nosFn, status, ...
                nSpikeUsed, nNoUsed, nMeta, nAgree}; %#ok<SAGROW>
            continue;
        end

        nMeta = height(Tm);
        idxAgree = find(round(Tm.agree) == 1);
        nAgree = numel(idxAgree);

        % In our reconstructed/event-export format:
        % spike_events.mat stores only accepted/valid spike epochs.
        nE = size(Xs,1);
        nUse = min(nE, nAgree);

        if nE ~= nAgree
            warning('Spike/meta mismatch for %s (%s): epochs=%d, meta agree==1=%d -> using %d', ...
                rname, base, nE, nAgree, nUse);
            status = 'SPIKE_META_MISMATCH';
        end

        if nUse == 0
            warning('No usable spike epochs for %s | %s', rname, base);
            Rows(end+1,:) = {rname, base, spkFn, metaFn, nosFn, status, ...
                nSpikeUsed, nNoUsed, nMeta, nAgree}; %#ok<SAGROW>
            continue;
        end

        Xs_use = Xs(1:nUse,:,:);
        samp   = round(Tm.sample(idxAgree(1:nUse)));
        code   = round(Tm.code(idxAgree(1:nUse)));
        agree  = ones(nUse,1);

        sid = extract_subject_id(base);

        %% ---------- Add spike candidates into shared consensus store -------
        for ii = 1:nUse

            key = make_key(base, samp(ii), code(ii));

            if ~isKey(key2idx, key)

                idx = numel(X_cell) + 1;
                key2idx(key) = idx;

                xi = squeeze(Xs_use(ii,:,:));  % C x T

                if opts.useSingle
                    xi = single(xi);
                end

                X_cell{idx,1}    = xi;
                subj_cell{idx,1} = sid;
                base_cell{idx,1} = base;
                sample_v(idx,1)  = samp(ii);
                code_v(idx,1)    = code(ii);
                Fs_v(idx,1)      = Fs_s;

                labelsMat(idx,:) = NaN(1, numel(opts.raters));
                labelsMat(idx,r) = agree(ii);

            else

                idx = key2idx(key);

                if ~isnan(labelsMat(idx,r))
                    warning('Duplicate label for same rater/key: %s | %s', rname, key);
                end

                labelsMat(idx,r) = agree(ii);
            end
        end

        nSpikeUsed = nUse;

        %% ---------- Add no-spike epochs as negatives ----------------------
        if opts.addNoSpikeEpochs && isfile(nosFn)

            try
                [En, metan] = load_event_mat_flexible(nosFn);
            catch ME
                warning('Could not load no-spike MAT: %s\n%s', nosFn, ME.message);
                Rows(end+1,:) = {rname, base, spkFn, metaFn, nosFn, 'BAD_NOSPIKE_MAT', ...
                    nSpikeUsed, nNoUsed, nMeta, nAgree}; %#ok<SAGROW>
                continue;
            end

            if ndims(En) ~= 3
                warning('No-spike events are not 3D: %s', nosFn);
            else

                cropIdx = opts.cropSamp(opts.cropSamp <= size(En,2));
                En = En(:, cropIdx, :);

                Xn = permute(En, [1 3 2]);  % N x channels x time

                if size(Xn,2) ~= opts.expectedNChan || size(Xn,3) ~= opts.expectedNTime
                    warning('Unexpected no-spike size for %s: Xn = [%d x %d x %d]', ...
                        base, size(Xn,1), size(Xn,2), size(Xn,3));
                else

                    if opts.useSingle
                        Xn = single(Xn);
                    end

                    Fs_n = getFs(metan, opts.Fs_default);
                    nsCenters = get_nospike_centers(metan, size(Xn,1));

                    keepNo = true(size(Xn,1),1);

                    if opts.deduplicateNoSpike
                        for jj = 1:size(Xn,1)

                            if ~isnan(nsCenters(jj))
                                nsKey = sprintf('%s|%d', base, round(nsCenters(jj)));
                            else
                                % If no center info exists, do not deduplicate this row.
                                nsKey = sprintf('%s|%s|row%d', base, rname, jj);
                            end

                            if isKey(nospikeSeen, nsKey)
                                keepNo(jj) = false;
                            else
                                nospikeSeen(nsKey) = true;
                            end
                        end
                    end

                    Xn = Xn(keepNo,:,:);
                    nNo = size(Xn,1);

                    if nNo > 0
                        X_ns_cell{end+1,1} = Xn; %#ok<SAGROW>
                        subj_ns = [subj_ns; repmat({sid}, nNo, 1)]; %#ok<AGROW>
                        base_ns = [base_ns; repmat({base}, nNo, 1)]; %#ok<AGROW>
                        Fs_ns   = [Fs_ns; repmat(Fs_n, nNo, 1)]; %#ok<AGROW>
                    end

                    nNoUsed = nNo;
                end
            end
        end

        Rows(end+1,:) = {rname, base, spkFn, metaFn, nosFn, status, ...
            nSpikeUsed, nNoUsed, nMeta, nAgree}; %#ok<SAGROW>
    end
end

%% ======================== FILE READ SUMMARY ==============================

if isempty(Rows)
    error('No files were read. Check root path and input folders.');
end

Tpairs = cell2table(Rows, 'VariableNames', ...
    {'Rater','Base','SpikeFile','MetaFile','NoSpikeFile','Status', ...
     'NSpikeUsed','NNoSpikeUsed','NMetaRows','NAgreeMeta'});

disp(Tpairs);

summaryCsv = fullfile(opts.root, 'Step4_file_read_summary.csv');
writetable(Tpairs, summaryCsv);
fprintf('\nSaved file-read summary:\n%s\n', summaryCsv);

fprintf('\nUnique spike candidates pooled across raters: %d\n', numel(X_cell));

if isempty(X_cell)
    error('No spike candidates were pooled. Check meta agree labels and spike MAT files.');
end

%% ======================== BUILD CONSENSUS LABELS =========================

nCand = numel(X_cell);
nReviewed = sum(~isnan(labelsMat), 2);

keep = nReviewed >= opts.minRatersForConsensus;

cons = NaN(nCand,1);
agreeCount = NaN(nCand,1);
rejectCount = NaN(nCand,1);

for i = 1:nCand

    labs = labelsMat(i, ~isnan(labelsMat(i,:)));

    if isempty(labs)
        continue;
    end

    agreeCount(i)  = sum(labs == 1);
    rejectCount(i) = sum(labs == 0);

    switch lower(opts.consensusRule)

        case 'majority'
            m = mean(labs);

            if m > 0.5
                cons(i) = 1;
            elseif m < 0.5
                cons(i) = 0;
            else
                cons(i) = NaN;  % tie -> drop
            end

        otherwise
            error('Unknown consensus rule: %s', opts.consensusRule);
    end
end

keep = keep & ~isnan(cons);

idxKeep = find(keep);

fprintf('\nConsensus candidates kept: %d / %d\n', numel(idxKeep), nCand);
fprintf('Consensus positives: %d\n', sum(cons(idxKeep)==1));
fprintf('Consensus negatives from candidate labels: %d\n', sum(cons(idxKeep)==0));

if isempty(idxKeep)
    error('No consensus candidates kept. Check minRatersForConsensus or meta labels.');
end

%% ======================== CONVERT CONSENSUS X ============================

% Nkeep = numel(idxKeep);
% C = size(X_cell{idxKeep(1)}, 1);
% T = size(X_cell{idxKeep(1)}, 2);
% 
% X_cons = zeros(Nkeep, C, T, 'like', X_cell{idxKeep(1)});
% 
% for ii = 1:Nkeep
%     X_cons(ii,:,:) = X_cell{idxKeep(ii)};
% end
% 
% y_cons = cons(idxKeep);
% 
% subj_cons = subj_cell(idxKeep);
% base_cons = base_cell(idxKeep);
% Fs_cons   = Fs_v(idxKeep);

%% ======================== LOW-MEMORY BUILD + SAVE ========================
% Avoid creating X_cons and allX in RAM.
% Saves large arrays as top-level variables:
%   X_train, X_val, X_test
% and saves metadata in dataset struct.

fprintf('\nUsing low-memory dataset writer...\n');

%% ---- Positive spike metadata ----
posN = numel(idxKeep);

C = size(X_cell{idxKeep(1)}, 1);
T = size(X_cell{idxKeep(1)}, 2);

y_pos    = cons(idxKeep);
subj_pos = subj_cell(idxKeep);
base_pos = base_cell(idxKeep);
Fs_pos   = Fs_v(idxKeep);

%% ---- Flatten no-spike cell chunks without concatenating data ----
nsChunk = [];
nsRow   = [];

if opts.addNoSpikeEpochs && exist('X_ns_cell','var') && ~isempty(X_ns_cell)
    for cc = 1:numel(X_ns_cell)
        nThis = size(X_ns_cell{cc}, 1);
        nsChunk = [nsChunk; repmat(cc, nThis, 1)]; %#ok<AGROW>
        nsRow   = [nsRow; (1:nThis)'];             %#ok<AGROW>
    end
end

nsN = numel(nsChunk);

if nsN > 0
    if numel(subj_ns) ~= nsN
        error('No-spike metadata mismatch: numel(subj_ns)=%d but nsN=%d', numel(subj_ns), nsN);
    end

    ally    = [y_pos(:); zeros(nsN,1)];
    subj_id = [subj_pos(:); subj_ns(:)];
    base_id = [base_pos(:); base_ns(:)];
    Fs_list = [Fs_pos(:); Fs_ns(:)];
else
    ally    = y_pos(:);
    subj_id = subj_pos(:);
    base_id = base_pos(:);
    Fs_list = Fs_pos(:);
end

fprintf('Positive spike epochs: %d\n', posN);
fprintf('No-spike epochs available: %d\n', nsN);
fprintf('Total epochs: %d\n', numel(ally));
fprintf('Class counts: spike=%d, nospike=%d\n', sum(ally==1), sum(ally==0));

if sum(ally==0) == 0
    warning('No no-spike epochs were found. Dataset will be spike-only.');
end

%% ---- Determine Fs ----
Fs = mode(round(Fs_list(~isnan(Fs_list))));

if isempty(Fs) || isnan(Fs)
    Fs = opts.Fs_default;
end

%% ---- Subject split before building arrays ----
[isTrain, isVal, isTest, trainSubj, valSubj, testSubj] = ...
    make_subject_split(subj_id, ally, opts.splitRatio, opts.rngSeed);

fprintf('\nSubject split:\n');
fprintf('Subjects: train=%d, val=%d, test=%d\n', ...
    numel(trainSubj), numel(valSubj), numel(testSubj));

fprintf('Trials:    train=%d, val=%d, test=%d\n', ...
    sum(isTrain), sum(isVal), sum(isTest));

fprintf('Train class counts: spike=%d, nospike=%d\n', ...
    sum(ally(isTrain)==1), sum(ally(isTrain)==0));
fprintf('Val class counts:   spike=%d, nospike=%d\n', ...
    sum(ally(isVal)==1), sum(ally(isVal)==0));
fprintf('Test class counts:  spike=%d, nospike=%d\n', ...
    sum(ally(isTest)==1), sum(ally(isTest)==0));

%% ---- Compute train normalization statistics without concatenating data ----
fprintf('\nComputing TRAIN normalization statistics in low-memory mode...\n');

trainGlobalIdx = find(isTrain);

sumC   = zeros(C,1);
sumSqC = zeros(C,1);
nObs   = 0;

for ii = 1:numel(trainGlobalIdx)

    g = trainGlobalIdx(ii);

    xi = fetch_epoch_from_store(g, posN, idxKeep, X_cell, X_ns_cell, nsChunk, nsRow);
    xi = double(xi);  % C x T

    sumC   = sumC   + sum(xi, 2);
    sumSqC = sumSqC + sum(xi.^2, 2);
    nObs   = nObs + size(xi,2);

    if mod(ii,1000) == 0
        fprintf('  stats: %d / %d train epochs\n', ii, numel(trainGlobalIdx));
    end
end

mu = sumC ./ nObs;

varC = (sumSqC - (sumC.^2 ./ nObs)) ./ max(nObs - 1, 1);
varC(varC < 0) = 0;

sig = sqrt(varC);
sig(sig == 0 | isnan(sig)) = 1;
mu(isnan(mu)) = 0;

mu_s  = single(mu);
sig_s = single(sig);

fprintf('Computed mu/sig for %d channels using %d observations per channel.\n', C, nObs);

%% ---- Prepare output file ----
if isfile(opts.outFile)
    backupFile = [opts.outFile '.bak_' datestr(now,'yyyymmdd_HHMMSS')];
    movefile(opts.outFile, backupFile);
    fprintf('Existing output backed up to:\n%s\n', backupFile);
end

nTrain = sum(isTrain);
nVal   = sum(isVal);
nTest  = sum(isTest);

y_train = ally(isTrain);
y_val   = ally(isVal);
y_test  = ally(isTest);

subject_id_train = subj_id(isTrain);
subject_id_val   = subj_id(isVal);
subject_id_test  = subj_id(isTest);

base_train = base_id(isTrain);
base_val   = base_id(isVal);
base_test  = base_id(isTest);

consensus = struct();
consensus.raters = opts.raters;
consensus.minRaters = opts.minRatersForConsensus;
consensus.rule = opts.consensusRule;
consensus.labelsMat = labelsMat(idxKeep,:);
consensus.nReviewed = nReviewed(idxKeep);
consensus.agreeCount = agreeCount(idxKeep);
consensus.rejectCount = rejectCount(idxKeep);
consensus.sample = sample_v(idxKeep);
consensus.code = code_v(idxKeep);

dataset = struct();
dataset.storage = 'low_memory_top_level_arrays';
dataset.note = 'Large arrays are saved as top-level variables: X_train, X_val, X_test.';
dataset.X_train_variable = 'X_train';
dataset.X_val_variable   = 'X_val';
dataset.X_test_variable  = 'X_test';

dataset.y_train = y_train;
dataset.y_val   = y_val;
dataset.y_test  = y_test;

dataset.subject_id_train = subject_id_train;
dataset.subject_id_val   = subject_id_val;
dataset.subject_id_test  = subject_id_test;

dataset.base_train = base_train;
dataset.base_val   = base_val;
dataset.base_test  = base_test;

dataset.trainSubj = trainSubj;
dataset.valSubj   = valSubj;
dataset.testSubj  = testSubj;

dataset.mu = mu;
dataset.sig = sig;
dataset.Fs = Fs;
dataset.consensus = consensus;
dataset.options = opts;
dataset.fileSummary = Tpairs;
dataset.buildDate = datestr(now);

Tsummary = table();

Tsummary.N_total = numel(ally);
Tsummary.N_spike = sum(ally==1);
Tsummary.N_nospike = sum(ally==0);

Tsummary.N_train = nTrain;
Tsummary.N_train_spike = sum(y_train==1);
Tsummary.N_train_nospike = sum(y_train==0);

Tsummary.N_val = nVal;
Tsummary.N_val_spike = sum(y_val==1);
Tsummary.N_val_nospike = sum(y_val==0);

Tsummary.N_test = nTest;
Tsummary.N_test_spike = sum(y_test==1);
Tsummary.N_test_nospike = sum(y_test==0);

Tsummary.N_subjects_total = numel(unique(subj_id));
Tsummary.N_subjects_train = numel(trainSubj);
Tsummary.N_subjects_val = numel(valSubj);
Tsummary.N_subjects_test = numel(testSubj);

Tsummary.N_channels = C;
Tsummary.N_time = T;
Tsummary.Fs = Fs;

save(opts.outFile, ...
    'dataset', ...
    'y_train','y_val','y_test', ...
    'subject_id_train','subject_id_val','subject_id_test', ...
    'base_train','base_val','base_test', ...
    'trainSubj','valSubj','testSubj', ...
    'mu','sig','Fs','consensus','Tsummary', ...
    '-v7.3');

fprintf('\nCreated metadata file:\n%s\n', opts.outFile);

%% ---- Create large arrays on disk without allocating them in RAM ----
m = matfile(opts.outFile, 'Writable', true);

m.X_train(nTrain, C, T) = single(0);
m.X_val(nVal, C, T)     = single(0);
m.X_test(nTest, C, T)   = single(0);

% fprintf('\nWriting normalized epochs directly to disk...\n');
% 
% iTrain = 0;
% iVal   = 0;
% iTest  = 0;
% 
% mu3  = reshape(mu_s,  [], 1);
% sig3 = reshape(sig_s, [], 1);
% 
% for g = 1:numel(ally)
% 
%     xi = fetch_epoch_from_store(g, posN, idxKeep, X_cell, X_ns_cell, nsChunk, nsRow);
%     xi = single(xi);   % C x T
% 
%     xi = (xi - mu3) ./ sig3;
%     xi = reshape(xi, [1 C T]);
% 
%     if isTrain(g)
%         iTrain = iTrain + 1;
%         m.X_train(iTrain,:,:) = xi;
% 
%     elseif isVal(g)
%         iVal = iVal + 1;
%         m.X_val(iVal,:,:) = xi;
% 
%     elseif isTest(g)
%         iTest = iTest + 1;
%         m.X_test(iTest,:,:) = xi;
%     end
% 
%     if mod(g,1000) == 0
%         fprintf('  wrote %d / %d epochs\n', g, numel(ally));
%     end
% end
% 
% fprintf('\nFinished writing arrays.\n');
% fprintf('Wrote train=%d, val=%d, test=%d\n', iTrain, iVal, iTest);

%% ---- Write large arrays to disk in batches ------------------------------
fprintf('\nWriting normalized epochs directly to disk in batches...\n');

batchSize = 256;   % try 256; if memory issue, use 128 or 64

bufTrain = zeros(batchSize, C, T, 'single');
bufVal   = zeros(batchSize, C, T, 'single');
bufTest  = zeros(batchSize, C, T, 'single');

nBufTrain = 0;
nBufVal   = 0;
nBufTest  = 0;

iTrain = 0;
iVal   = 0;
iTest  = 0;

writeTrainStart = 1;
writeValStart   = 1;
writeTestStart  = 1;

% mu_s and sig_s are C x 1
for g = 1:numel(ally)

    xi = fetch_epoch_from_store(g, posN, idxKeep, X_cell, X_ns_cell, nsChunk, nsRow);
    xi = single(xi);   % C x T

    % Normalize: C x T
    xi = bsxfun(@rdivide, bsxfun(@minus, xi, mu_s), sig_s);

    % Store as 1 x C x T
    xi = reshape(xi, [1 C T]);

    if isTrain(g)

        nBufTrain = nBufTrain + 1;
        bufTrain(nBufTrain,:,:) = xi;

        if nBufTrain == batchSize
            m.X_train(writeTrainStart:writeTrainStart+nBufTrain-1,:,:) = bufTrain;
            writeTrainStart = writeTrainStart + nBufTrain;
            iTrain = iTrain + nBufTrain;
            nBufTrain = 0;
        end

    elseif isVal(g)

        nBufVal = nBufVal + 1;
        bufVal(nBufVal,:,:) = xi;

        if nBufVal == batchSize
            m.X_val(writeValStart:writeValStart+nBufVal-1,:,:) = bufVal;
            writeValStart = writeValStart + nBufVal;
            iVal = iVal + nBufVal;
            nBufVal = 0;
        end

    elseif isTest(g)

        nBufTest = nBufTest + 1;
        bufTest(nBufTest,:,:) = xi;

        if nBufTest == batchSize
            m.X_test(writeTestStart:writeTestStart+nBufTest-1,:,:) = bufTest;
            writeTestStart = writeTestStart + nBufTest;
            iTest = iTest + nBufTest;
            nBufTest = 0;
        end
    end

    if mod(g,1000) == 0
        fprintf('  prepared/wrote %d / %d epochs\n', g, numel(ally));
    end
end

%% ---- Flush remaining partial batches -----------------------------------

if nBufTrain > 0
    m.X_train(writeTrainStart:writeTrainStart+nBufTrain-1,:,:) = bufTrain(1:nBufTrain,:,:);
    iTrain = iTrain + nBufTrain;
end

if nBufVal > 0
    m.X_val(writeValStart:writeValStart+nBufVal-1,:,:) = bufVal(1:nBufVal,:,:);
    iVal = iVal + nBufVal;
end

if nBufTest > 0
    m.X_test(writeTestStart:writeTestStart+nBufTest-1,:,:) = bufTest(1:nBufTest,:,:);
    iTest = iTest + nBufTest;
end

clear bufTrain bufVal bufTest xi

fprintf('\nFinished writing arrays in batches.\n');
fprintf('Wrote train=%d, val=%d, test=%d\n', iTrain, iVal, iTest);

%% ---- Save dataset summary CSV ----
summaryFile = fullfile(opts.root, 'Step4_dataset_summary.csv');
writetable(Tsummary, summaryFile);
fprintf('Saved dataset summary:\n%s\n', summaryFile);

%% ---- Quick file check ----
fprintf('\nVariables in saved file:\n');
whos('-file', opts.outFile)

m = matfile(opts.outFile);

fprintf('\nSaved array sizes:\n');
fprintf('X_train: [%s]\n', num2str(size(m, 'X_train')));
fprintf('X_val:   [%s]\n', num2str(size(m, 'X_val')));
fprintf('X_test:  [%s]\n', num2str(size(m, 'X_test')));

fprintf('\nStep 4 complete in low-memory mode.\n');

%% ========================================================================
% Local helper functions
% ========================================================================

function [E, M] = load_event_mat_flexible(fn)

    S = load(fn);

    if isfield(S, 'events')
        E = S.events;
    elseif isfield(S, 'Es')
        E = S.Es;
    elseif isfield(S, 'X')
        E = S.X;
    else
        f = fieldnames(S);
        E = [];

        for i = 1:numel(f)
            v = S.(f{i});
            if isnumeric(v) && ndims(v) == 3
                E = v;
                break;
            end
        end

        if isempty(E)
            error('No 3D numeric event variable found in %s', fn);
        end
    end

    if isfield(S, 'metas')
        M = S.metas;
    elseif isfield(S, 'meta')
        M = S.meta;
    else
        M = struct();
    end
end

function T = read_meta_flexible(fn)

    try
        imp = detectImportOptions(fn, 'FileType', 'text');
        T = readtable(fn, imp);
    catch
        T = readtable(fn, 'FileType', 'text', 'Delimiter', '\t');
    end

    if isempty(T)
        return;
    end

    oldNames = T.Properties.VariableNames;
    normNames = lower(regexprep(oldNames, '[^a-zA-Z0-9]', ''));

    for i = 1:numel(normNames)

        switch normNames{i}

            case {'sample','sampleidx','sampleindex','eventsample','latency'}
                T.Properties.VariableNames{i} = 'sample';

            case {'agree','accepted','accept','include','included','isspike','spike','label','decision'}
                T.Properties.VariableNames{i} = 'agree';

            case {'code','eventcode','eventid','trigger','value'}
                T.Properties.VariableNames{i} = 'code';

            case {'sampleds','sampledownsampled','sample200hz'}
                T.Properties.VariableNames{i} = 'sample_ds';
        end
    end

    if ismember('sample', T.Properties.VariableNames)
        T.sample = col_to_double(T.sample);
    end

    if ismember('agree', T.Properties.VariableNames)
        T.agree = col_to_agree(T.agree);
    end

    if ismember('code', T.Properties.VariableNames)
        T.code = col_to_double(T.code);
    end
end

function x = col_to_double(v)

    if isnumeric(v)
        x = double(v);
    elseif iscell(v)
        x = str2double(string(v));
    elseif isstring(v) || ischar(v) || iscategorical(v)
        x = str2double(string(v));
    else
        x = double(v);
    end

    x = x(:);
end

function x = col_to_agree(v)

    if isnumeric(v)
        x = double(v(:));
        return;
    end

    s = lower(strtrim(string(v)));
    x = nan(numel(s),1);

    yesVals = ["1","yes","y","true","t","spike","agree","accepted","accept","include","included"];
    noVals  = ["0","no","n","false","f","nospike","reject","rejected","exclude","excluded"];

    x(ismember(s, yesVals)) = 1;
    x(ismember(s, noVals)) = 0;

    nums = str2double(s);
    useNum = isnan(x) & ~isnan(nums);
    x(useNum) = nums(useNum);
end

function Fs = getFs(meta, Fs_default)

    Fs = Fs_default;

    if isstruct(meta) && isfield(meta, 'Fs') && ~isempty(meta.Fs)
        Fs = double(meta.Fs);
        Fs = Fs(1);
    end
end

function sid = extract_subject_id(base)

    tok = regexp(base, '^(.*)_Run\d+', 'tokens', 'once');

    if ~isempty(tok)
        sid = tok{1};
    else
        sid = base;
    end
end

function key = make_key(base, sample, code)

    if isnan(sample)
        sample = -1;
    end

    if isnan(code)
        code = -1;
    end

    key = sprintf('%s|%d|%d', base, round(sample), round(code));
end

function centers = get_nospike_centers(meta, n)

    centers = nan(n,1);

    if isstruct(meta)

        if isfield(meta, 'sample_ds') && numel(meta.sample_ds) >= n
            centers = double(meta.sample_ds(1:n));
            centers = centers(:);
            return;
        end

        if isfield(meta, 'sample') && numel(meta.sample) >= n
            centers = double(meta.sample(1:n));
            centers = centers(:);
            return;
        end
    end
end

function [isTrain, isVal, isTest, trainSubj, valSubj, testSubj] = ...
    make_subject_split(subj_id, y, splitRatio, seed)

    rng(seed);

    subjects = unique(subj_id, 'stable');
    nSubj = numel(subjects);

    if nSubj < 3
        error('Need at least 3 subjects for train/val/test split. Found %d.', nSubj);
    end

    nTrain = round(splitRatio(1) * nSubj);
    nVal   = round(splitRatio(2) * nSubj);

    nTrain = max(1, nTrain);
    nVal   = max(1, nVal);

    if nTrain + nVal >= nSubj
        nVal = max(1, nSubj - nTrain - 1);
    end

    maxTry = 1000;
    best = [];

    for tr = 1:maxTry

        idx = randperm(nSubj);

        trainSubj = subjects(idx(1:nTrain));
        valSubj   = subjects(idx(nTrain+1:nTrain+nVal));
        testSubj  = subjects(idx(nTrain+nVal+1:end));

        isTrain = ismember(subj_id, trainSubj);
        isVal   = ismember(subj_id, valSubj);
        isTest  = ismember(subj_id, testSubj);

        okTrain = has_both_classes(y(isTrain));
        okVal   = has_both_classes(y(isVal));
        okTest  = has_both_classes(y(isTest));

        score = double(okTrain) + double(okVal) + double(okTest);

        if isempty(best) || score > best.score
            best.trainSubj = trainSubj;
            best.valSubj = valSubj;
            best.testSubj = testSubj;
            best.isTrain = isTrain;
            best.isVal = isVal;
            best.isTest = isTest;
            best.score = score;
        end

        if okTrain && okVal && okTest
            return;
        end
    end

    warning('Could not find a split where train/val/test all contain both classes. Using best available split.');

    trainSubj = best.trainSubj;
    valSubj = best.valSubj;
    testSubj = best.testSubj;

    isTrain = best.isTrain;
    isVal = best.isVal;
    isTest = best.isTest;
end

function tf = has_both_classes(y)

    tf = any(y==1) && any(y==0);
end

function [Xn, mu, sig] = normalize_by_train(X)

    [~, C, ~] = size(X);

    tmp = permute(X, [2 1 3]);  % C x N x T
    tmp = reshape(tmp, C, []);  % C x (N*T)

    mu = mean(tmp, 2, 'omitnan');
    sig = std(tmp, 0, 2, 'omitnan');

    sig(sig == 0 | isnan(sig)) = 1;
    mu(isnan(mu)) = 0;

    Xn = apply_norm(X, mu, sig);
end

function Xn = apply_norm(X, mu, sig)

    mu3 = reshape(mu, 1, [], 1);
    sig3 = reshape(sig, 1, [], 1);

    Xn = (X - mu3) ./ sig3;
end

function xi = fetch_epoch_from_store(g, posN, idxKeep, X_cell, X_ns_cell, nsChunk, nsRow)

    if g <= posN
        xi = X_cell{idxKeep(g)};   % C x T
    else
        j = g - posN;
        cc = nsChunk(j);
        rr = nsRow(j);

        xi = squeeze(X_ns_cell{cc}(rr,:,:));  % C x T
    end

    if ndims(xi) ~= 2
        error('Fetched epoch is not 2D. Size = [%s]', num2str(size(xi)));
    end
end