%% ------------------------------------------------------------------------
% Build Spike/NoSpike dataset with multi-rater consensus from meta .txt
% + leak-free split + normalization
% -------------------------------------------------------------------------
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/function')
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/MEGneto/func')

root   = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files';
raters = {'Adi','Josh','Manoj','Pradeep'};

% -------------------- User settings --------------------------------------
cropSamp    = 1:306;     % keep first 306 sensors (your current choice)
useSingle   = true;      % store X as single to save RAM
Fs_default  = 200;       % fallback if Fs not stored
rng(1);                  % reproducible split

minRatersForConsensus = 2;       % keep candidates reviewed by >= this many raters
consensusRule         = 'majority'; % 'majority' (ties -> drop)
addNoSpikeEpochs      = true;    % also include *_nospike_events.mat as extra negatives

% -------------------- Candidate store (keyed by base+sample+code) --------
key2idx = containers.Map('KeyType','char','ValueType','double');

X_cell   = {};           % each cell: 1xCxT
subj_cell= {};
base_cell= {};
sample_v = [];           % Nx1
code_v   = [];           % Nx1
Fs_v     = [];           % Nx1 (best effort)
labelsMat = NaN(0, numel(raters));  % NxR labels per rater (0/1)

% For bookkeeping of nospike epochs (added directly later)
X_ns = [];
y_ns = [];
subj_ns = {};
base_ns = {};
Fs_ns = [];

pairs = {};  % summary table per rater/run

for r = 1:numel(raters)
    rname = raters{r};
    d = fullfile(root, rname);
    
    spk = dir(fullfile(d, '*_spike_events.mat'));
    
    for k = 1:numel(spk)
        spk_fn = fullfile(spk(k).folder, spk(k).name);
        base   = erase(spk(k).name, '_spike_events.mat');
        
        nos_fn = fullfile(d, [base '_nospike_events.mat']);
        meta_fn = fullfile(d, 'meta_files', sprintf('%s_%s_meta_data.txt', base, rname));
        
        % ----- Load spike epochs -----
        [Es, metas] = load_epochs_mat(spk_fn);   % Es: trials x time x chan (expected)
        Fs_s = getFs(metas, Fs_default);
        
        % Crop time safely
        Es = Es(:, cropSamp(cropSamp <= size(Es,2)), :);
        
        % Convert to N×C×T
        Xs = permute(Es, [1 3 2]);
        
        % ----- Read meta labels (agree 0/1) -----
        meta = read_meta_txt(meta_fn, Fs_s, false);
        Tm   = meta.table;
        
        if isempty(Tm) || ~all(ismember({'sample','agree','code'}, Tm.Properties.VariableNames))
            fprintf('Skipping (no usable meta): %s\n', meta_fn);
            continue;
        end
        
        %         % Align meta rows to epochs (assumes same order; trims if mismatch)
        %         nE = size(Xs,1);
        %         nM = height(Tm);
        %         nUse = min(nE, nM);
        %         if nE ~= nM
        %             warning('Epoch/meta mismatch for %s (%s): epochs=%d, meta=%d -> using %d', ...
        %                 rname, base, nE, nM, nUse);
        %         end
        %
        %         Xs_use = Xs(1:nUse,:,:);
        %         samp   = round(Tm.sample(1:nUse));
        %         agree  = round(Tm.agree(1:nUse));
        %         code   = round(Tm.code(1:nUse));
        
        % If spike_events.mat stores APPROVED spikes only:
        idx1 = find(round(Tm.agree)==1);
        n1 = numel(idx1);
        nE = size(Xs,1);
        
        nUse = min(nE, n1);
        if nE ~= n1
            warning('Spike-only export mismatch for %s (%s): epochs=%d but meta agree==1 count=%d -> using %d', ...
                rname, base, nE, n1, nUse);
        end
        
        Xs_use = Xs(1:nUse,:,:);
        samp   = round(Tm.sample(idx1(1:nUse)));
        agree  = ones(nUse,1);              % since these are the approved spikes
        code   = round(Tm.code(idx1(1:nUse)));
        
        
        sid = extract_subject_id(base);
        
        % ----- Add each candidate (base+sample+code) into shared store -----
        for i = 1:nUse
            key = make_key(base, samp(i), code(i));
            
            if ~isKey(key2idx, key)
                % new candidate
                idx = numel(X_cell) + 1;
                key2idx(key) = idx;
                
                xi = squeeze(Xs_use(i,:,:));      % C×T
                if useSingle, xi = single(xi); end
                
                X_cell{idx,1}    = xi;
                subj_cell{idx,1} = sid;
                base_cell{idx,1} = base;
                sample_v(idx,1)  = samp(i);
                code_v(idx,1)    = code(i);
                Fs_v(idx,1)      = Fs_s;
                
                labelsMat(idx, :) = NaN(1, numel(raters));
                labelsMat(idx, r) = agree(i);
            else
                % existing candidate, just fill this rater's label
                idx = key2idx(key);
                labelsMat(idx, r) = agree(i);
                
                % Optional sanity check: same epoch data?
                % (commented; can be expensive)
                % xi = squeeze(Xs_use(i,:,:));
                % if norm(double(xi(:)) - double(X_cell{idx}(:))) > 1e-6
                %     warning('Data mismatch for same key: %s (%s)', key, rname);
                % end
            end
        end
        
        % ----- Optionally add nospike epochs as extra negatives (no consensus needed) -----
        nNo = 0;
        if addNoSpikeEpochs && isfile(nos_fn)
            [En, metan] = load_epochs_mat(nos_fn);
            Fs_n = getFs(metan, Fs_default);
            
            En = En(:, cropSamp(cropSamp <= size(En,2)), :);
            Xn = permute(En, [1 3 2]);  % N×C×T
            if useSingle, Xn = single(Xn); end
            
            nNo = size(Xn,1);
            X_ns   = cat(1, X_ns, Xn);
            y_ns   = cat(1, y_ns, zeros(nNo,1));
            subj_ns= [subj_ns; repmat({sid}, nNo, 1)];
            base_ns= [base_ns; repmat({base}, nNo, 1)];
            Fs_ns  = [Fs_ns;  repmat(Fs_n, nNo, 1)]; %#ok<AGROW>
        end
        
        pairs(end+1,:) = {rname, base, spk_fn, meta_fn, nUse, nNo, size(Xs,2), size(Xs,3)}; 
    end
end

Tpairs = cell2table(pairs, 'VariableNames', ...
    {'Rater','Base','SpikeFile','MetaFile','nCandidatesUsed','nNoSpikeAdded','nChan','nTime'});
disp(Tpairs)

fprintf('Unique candidates pooled across raters: %d\n', numel(X_cell));

%% -------------------- Build consensus-labeled dataset --------------------
nCand = numel(X_cell);
nReviewed = sum(~isnan(labelsMat), 2);

keep = nReviewed >= minRatersForConsensus;

cons = NaN(nCand,1);
agreeCount = NaN(nCand,1);
rejectCount = NaN(nCand,1);

for i = 1:nCand
    labs = labelsMat(i, ~isnan(labelsMat(i,:)));
    if isempty(labs), continue; end
    
    agreeCount(i)  = sum(labs==1);
    rejectCount(i) = sum(labs==0);
    
    switch lower(consensusRule)
        case 'majority'
            m = mean(labs);
            if m > 0.5
                cons(i) = 1;
            elseif m < 0.5
                cons(i) = 0;
            else
                cons(i) = NaN; % tie -> drop
            end
        otherwise
            error('Unknown consensusRule: %s', consensusRule);
    end
end

keep = keep & ~isnan(cons);

% Convert X_cell -> numeric 3D array (N×C×T)
idxKeep = find(keep);
Nkeep = numel(idxKeep);
C = size(X_cell{idxKeep(1)},1);
T = size(X_cell{idxKeep(1)},2);

X_cons = zeros(Nkeep, C, T, 'like', X_cell{idxKeep(1)});
for ii = 1:Nkeep
    X_cons(ii,:,:) = X_cell{idxKeep(ii)};
end

y_cons = cons(idxKeep);

subj_cons = subj_cell(idxKeep);
base_cons = base_cell(idxKeep);

% Combine with nospike extras (if enabled)
allX = X_cons;
ally = y_cons;
subj_id = subj_cons;
base_id = base_cons;
Fs_list = Fs_v(idxKeep);

if addNoSpikeEpochs && ~isempty(X_ns)
    allX   = cat(1, allX, X_ns);
    ally   = cat(1, ally, y_ns);
    subj_id= [subj_id; subj_ns];
    base_id= [base_id; base_ns];
    Fs_list= [Fs_list; Fs_ns]; 
end

fprintf('Final pooled dataset: X=[%d×%d×%d], y=[%d×1]\n', size(allX,1), size(allX,2), size(allX,3), numel(ally));
fprintf('Positives (Spike=1): %d, Negatives (0): %d\n', sum(ally==1), sum(ally==0));

%% -------------------- Quick sanity plots --------------------------------
Fs = mode(Fs_list(~isnan(Fs_list)));
if isempty(Fs), Fs = Fs_default; end
t = (0:size(allX,3)-1)/Fs;

% trial = 1;
% figure; plot(t, squeeze(allX(trial,1,:)));
% xlabel('Time (s)'); ylabel('Amplitude');
% title(sprintf('Trial %d, Channel 1, y=%d', trial, ally(trial)));
% grid on;

figure; imagesc(t, 1:size(allX,2), squeeze(allX(trial,:,:)));
axis tight; xlabel('Time (s)'); ylabel('Channel');
title('Trial heatmap (channels × time)'); colorbar;

%% -------------------- Leak-free split by SUBJECT ------------------------
subjects = unique(subj_id);
nSubj = numel(subjects);

idx = randperm(nSubj);
nTrain = round(0.70*nSubj);
nVal   = round(0.15*nSubj);

trainSubj = subjects(idx(1:nTrain));
valSubj   = subjects(idx(nTrain+1:nTrain+nVal));
testSubj  = subjects(idx(nTrain+nVal+1:end));

isTrain = ismember(subj_id, trainSubj);
isVal   = ismember(subj_id, valSubj);
isTest  = ismember(subj_id, testSubj);

fprintf('Subjects: train=%d, val=%d, test=%d\n', numel(trainSubj), numel(valSubj), numel(testSubj));
fprintf('Trials:    train=%d, val=%d, test=%d\n', sum(isTrain), sum(isVal), sum(isTest));

%% -------------------- Normalize using TRAIN stats only -------------------
[X_train_norm, mu, sig] = normalize_by_train(allX(isTrain,:,:));
X_val_norm  = apply_norm(allX(isVal,:,:),  mu, sig);
X_test_norm = apply_norm(allX(isTest,:,:), mu, sig);

y_train = ally(isTrain);
y_val   = ally(isVal);
y_test  = ally(isTest);

%% -------------------- Save DL-ready dataset ------------------------------
dataset = struct();
dataset.X_train = X_train_norm; dataset.y_train = y_train;
dataset.X_val   = X_val_norm;   dataset.y_val   = y_val;
dataset.X_test  = X_test_norm;  dataset.y_test  = y_test;

dataset.subject_id_train = subj_id(isTrain);
dataset.subject_id_val   = subj_id(isVal);
dataset.subject_id_test  = subj_id(isTest);

dataset.base_train = base_id(isTrain);
dataset.base_val   = base_id(isVal);
dataset.base_test  = base_id(isTest);

% keep multi-rater info for the consensus candidates (before adding nospike)
dataset.consensus = struct();
dataset.consensus.raters = raters;
dataset.consensus.minRaters = minRatersForConsensus;
dataset.consensus.rule = consensusRule;
dataset.consensus.labelsMat = labelsMat(idxKeep,:);         % Nkeep × R
dataset.consensus.nReviewed = nReviewed(idxKeep);
dataset.consensus.agreeCount = agreeCount(idxKeep);
dataset.consensus.rejectCount = rejectCount(idxKeep);
dataset.consensus.sample = sample_v(idxKeep);
dataset.consensus.code   = code_v(idxKeep);

dataset.mu = mu; dataset.sig = sig;
dataset.Fs = Fs;

save('dataset_DL_ready.mat', 'dataset', '-v7.3');
disp('Saved: dataset_DL_ready.mat');

%% -------------------- Local helper functions ----------------------------
function Fs = getFs(meta, Fs_default)
Fs = Fs_default;
if isstruct(meta) && isfield(meta,'Fs') && ~isempty(meta.Fs)
    Fs = double(meta.Fs);
    if numel(Fs) > 1, Fs = Fs(1); end
end
end

function sid = extract_subject_id(base)
tok = regexp(base, '^(.*)_Run\d+', 'tokens', 'once');
if ~isempty(tok), sid = tok{1};
else, sid = base;
end
end

function key = make_key(base, sample, code)
key = sprintf('%s|%d|%d', base, sample, code);
end

function [Xn, mu, sig] = normalize_by_train(X)
[~,C,~] = size(X);
tmp = permute(X, [2 1 3]);     % C×N×T
tmp = reshape(tmp, C, []);     % C×(N*T)
mu  = mean(tmp, 2, 'omitnan');
sig = std(tmp, 0, 2, 'omitnan');
sig(sig==0) = 1;
Xn = apply_norm(X, mu, sig);
end

function Xn = apply_norm(X, mu, sig)
mu3  = reshape(mu,  1, [], 1);
sig3 = reshape(sig, 1, [], 1);
Xn = (X - mu3) ./ sig3;
end


% %% -------------------- Subject-level then group-level pair summaries -----
% [RunPairStats, PairGroupRuns] = summarize_pairs_by_run(labelsMat, base_cell, raters);
% 
% disp('=== Run-level stats (pair x base/run) ===');
% disp(RunPairStats(1:min(30,height(RunPairStats)),:));
% 
% disp('=== Group-level stats aggregated across runs (per pair) ===');
% disp(PairGroupRuns);
% 
% writetable(RunPairStats,  'PairRunLevel.csv');
% writetable(PairGroupRuns, 'PairGroupLevel_byRun.csv');
% 
% 
% %% -------------------- Pairwise reviewer summary (ALL pairs) -------------
% doPairSummary = true;
% 
% if doPairSummary
%     [PairSummary, WorstRuns] = summarize_all_rater_pairs(labelsMat, base_cell, raters);
%     
%     disp('=== Pairwise summary across ALL candidates ===');
%     disp(PairSummary);
%     
%     disp('=== Worst disagreement runs per pair (top 10 each) ===');
%     disp(WorstRuns);
% end
% 
% doPairSummary = true;
% 
% if doPairSummary
%     [PairSummary, WorstRuns, FullRuns] = summarize_all_rater_pairs(labelsMat, base_cell, raters);
%     
%     disp('=== Full per-run list (all bases x all rater pairs) ===');
%     disp(FullRuns(1:min(30,height(FullRuns)),:));   % preview first 30
%     
%     %     disp('=== Worst disagreement runs per pair (top 10 each) ===');
%     %     disp(WorstRuns);
% end
% 
% %%
% for r = 1:numel(raters)
%     idx = ~all(isnan(labelsMat(:,r)),2);  % candidates with label from r
%     fprintf('%s: candidates labeled=%d, unique bases=%d\n', ...
%         raters{r}, sum(~isnan(labelsMat(:,r))), numel(unique(base_cell(idx))));
% end
% 
% for i=1:numel(raters)
%     for j=i+1:numel(raters)
%         a = labelsMat(:,i); b = labelsMat(:,j);
%         nO = sum(~isnan(a) & ~isnan(b));
%         fprintf('%s vs %s overlap candidates = %d\n', raters{i}, raters{j}, nO);
%     end
% end
% 
% for i=1:numel(raters)
%     for j=i+1:numel(raters)
%         bases = unique(base_cell);
%         nBaseOverlap = 0;
%         for b=1:numel(bases)
%             idxB = strcmp(base_cell, bases{b});
%             a = labelsMat(idxB,i); bb = labelsMat(idxB,j);
%             if any(~isnan(a) & ~isnan(bb))
%                 nBaseOverlap = nBaseOverlap + 1;
%             end
%         end
%         fprintf('%s vs %s overlap bases = %d\n', raters{i}, raters{j}, nBaseOverlap);
%     end
% end
% 
% % tolSamp = round(0.010 * Fs_s);   % 10 ms tolerance (adjust)
% % samp_key = tolSamp * round(samp(i)/tolSamp);  % quantize
% % key = make_key(base, samp_key, code(i));
% 
% %% -------------------- Individual-level two-reviewer plots (optional) ----
% doPlotTwoRater = true;      % toggle
% plotR1 = 'Manoj';
% plotR2 = 'Josh';
% 
% if doPlotTwoRater
%     i1 = find(strcmp(raters, plotR1));
%     i2 = find(strcmp(raters, plotR2));
%     if isempty(i1) || isempty(i2)
%         warning('plotR1/plotR2 not found in raters list.');
%     else
%         % Choose a base that has overlap for both reviewers
%         ubase = unique(base_cell);
%         bestBase = '';
%         bestN = 0;
%         
%         for bb = 1:numel(ubase)
%             b = ubase{bb};
%             idxB = find(strcmp(base_cell, b));
%             a1 = labelsMat(idxB, i1);
%             a2 = labelsMat(idxB, i2);
%             nOverlap = sum(~isnan(a1) & ~isnan(a2));
%             if nOverlap > bestN
%                 bestN = nOverlap;
%                 bestBase = b;
%             end
%         end
%         
%         if bestN == 0
%             warning('No overlapping reviews found for %s vs %s.', plotR1, plotR2);
%         else
%             fprintf('Plotting base with most overlap (%s vs %s): %s (overlap=%d)\n', ...
%                 plotR1, plotR2, bestBase, bestN);
%             
%             % Build per-base “Tm-like” table from your pooled store
%             idxB = find(strcmp(base_cell, bestBase));
%             Tb = table(sample_v(idxB), code_v(idxB), ...
%                 labelsMat(idxB,i1), labelsMat(idxB,i2), ...
%                 'VariableNames', {'sample','code',['agree_' plotR1],['agree_' plotR2]});
%             
%             % Add derived categories
%             a1 = Tb.(['agree_' plotR1]); a2 = Tb.(['agree_' plotR2]);
%             Tb.bothReviewed = ~isnan(a1) & ~isnan(a2);
%             Tb.bothAgree    = Tb.bothReviewed & (a1==1) & (a2==1);
%             Tb.bothReject   = Tb.bothReviewed & (a1==0) & (a2==0);
%             Tb.disagree     = Tb.bothReviewed & (a1~=a2);
%             
%             % Add time if Fs is known
%             Fs_here = mode(Fs_v(idxB));
%             if ~isempty(Fs_here) && ~isnan(Fs_here) && Fs_here > 0
%                 Tb.time_sec = double(Tb.sample) ./ double(Fs_here);
%             end
%             
%             % ----- Timeline plot -----
%             plot_two_rater_timeline(Tb, plotR1, plotR2, bestBase);
%             
%             % ----- Data category plots (mean templates + example traces) -----
%             % Pull X for this base from X_cell:
%             Xb = cell2mat(reshape(X_cell(idxB), 1, 1, [])); %#ok<NASGU>
%             % cell2mat won’t work directly for CxT cells; so do manual:
%             Nb = numel(idxB);
%             Cb = size(X_cell{idxB(1)},1);
%             TbT = size(X_cell{idxB(1)},2);
%             X3 = zeros(Nb, Cb, TbT, 'like', X_cell{idxB(1)});
%             for ii = 1:Nb
%                 X3(ii,:,:) = X_cell{idxB(ii)};
%             end
%             
%             plot_two_rater_data_categories_fromX(X3, Tb, plotR1, plotR2, Fs_here);
%         end
%     end
% end

