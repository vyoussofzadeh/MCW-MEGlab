%% ------------------------------------------------------------------------
% Step6_Batch_Run_Model_On_All_Continuous_MEG.m
%
% Batch apply trained spike/no-spike model to all continuous MEG FIF files.
%
% For each FIF:
%   1) Load continuous MEG
%   2) Sliding-window model prediction
%   3) Keep strong probability peaks only
%   4) Refine event time to local GFP peak
%   5) Deduplicate nearby events
%   6) Save:
%        *_model_candidates_gt75_refined_timepoints.csv
%        *_model_candidates_gt75_refined_timepoints.xlsx
%        *_model_candidates_gt75_refined_spk_DS.eve
%        *_model_candidates_gt75_refined_full.mat
%
% Excel/CSV contains only:
%   sample, time_sec
%
% .eve file follows old mbrowse-compatible format:
%   first_samp    first_samp/fs_hdr    0    0    test
%   event_sample  event_time_sec       0    5555
% -------------------------------------------------------------------------

clc; clear;

%% ======================== USER SETTINGS =================================

opts.ftRoot  = '/home/vyoussofzadeh/Desktop/research_workspace/tools/fieldtrip/fieldtrip_2022/';
opts.mneRoot = '/MEG_data/MEG_Tools/mne';

% Trained model from Step 5
% opts.modelFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Squiggles/Step5_lowmem_CNN1D_20260527_153010.mat';
% opts.modelFile = '/home/vyoussofzadeh/Data/DL_model/Step5_lowmem_CNN1D_20260602_181546.mat';
opts.modelFile = '/home/vyoussofzadeh/Data/DL_model/Step5_lowmem_CNN1D_bp5_50_m250_p500_20260604_134703.mat';


% If model file points to a dataset path that does not exist on this machine,
% set this to the correct dataset used for training. Otherwise leave empty.
opts.datasetFileOverride = '';

% Folder with FIF files
opts.rawDir = '/home/vyoussofzadeh/Data/DL_Spike/';

% File pattern
opts.filePattern = '*_raw_DS.fif';

% Search subfolders too?
opts.recursive = false;

% Output folder
% opts.outDir = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data/model_outputs_gt75';
opts.outDir = '/home/vyoussofzadeh/Data/DL_model_output//Data/model_outputs_gt75_bp5_50_m250_p500';

% Skip files if output MAT already exists
opts.skipExisting = true;

% Sliding-window settings
opts.stepSec = 0.050;                 % 50 ms step. Use 0.010 for finer timing.
% opts.modelWinSec = [-0.500 1.025];    % match dataset used for model training
opts.modelWinSec = [-0.250 0.500];

% Strong candidate selection
opts.probThreshold = 0.75;
opts.probSmoothSec = 0.050;
opts.minEventSepSec = 0.300;

% Local peak refinement
opts.refineToGfpPeak = true;
opts.refineWinSec = 0.150;
opts.refineDedupSec = 0.200;

% Event code for mbrowse/mne_browse_raw
opts.eventCode = 5555;

% Optional: process only part of each recording for testing
opts.useTimeRange = false;
opts.timeRangeSec = [0 300];

% Save complete probability trace in MAT file
opts.saveFullProbabilityTrace = true;

% Prediction mini-batch size
opts.miniBatchWindows = 128;

% Bandpass settings: must match BP model training
opts.useBandpass = true;
opts.bpFreq = [5 50];
opts.bpOrder = 4;

%% ======================== INITIALIZE ====================================

if ~exist(opts.outDir, 'dir')
    mkdir(opts.outDir);
end

if exist(opts.ftRoot, 'dir')
    addpath(opts.ftRoot);
    ft_defaults;
else
    error('FieldTrip folder not found: %s', opts.ftRoot);
end

if exist(opts.mneRoot, 'dir')
    addpath(genpath(opts.mneRoot));
end

rehash toolboxcache;

fprintf('\nChecking readers:\n');
which ft_preprocessing -all
which fiff_find_evoked -all

%% ======================== COLLECT FIF FILES ==============================

if opts.recursive
    L = dir(fullfile(opts.rawDir, '**', opts.filePattern));
else
    L = dir(fullfile(opts.rawDir, opts.filePattern));
end

if isempty(L)
    error('No FIF files found: %s', fullfile(opts.rawDir, opts.filePattern));
end

rawFiles = fullfile({L.folder}, {L.name})';

fprintf('\nFound %d FIF files:\n', numel(rawFiles));
for i = 1:numel(rawFiles)
    fprintf('[%03d] %s\n', i, rawFiles{i});
end

%% ======================== LOAD MODEL ONCE ================================

fprintf('\nLoading model:\n%s\n', opts.modelFile);

M = load(opts.modelFile);

if isfield(M, 'bestNet')
    net = M.bestNet;
elseif isfield(M, 'net')
    net = M.net;
elseif isfield(M, 'netBest')
    net = M.netBest;
else
    error('Model file does not contain bestNet, net, or netBest.');
end

if isfield(M, 'datasetFile')
    datasetFile = M.datasetFile;
else
    datasetFile = '';
end

if ~isempty(opts.datasetFileOverride)
    datasetFile = opts.datasetFileOverride;
end

if isempty(datasetFile) || ~isfile(datasetFile)
    error(['Could not find datasetFile used for normalization.\n' ...
           'Set opts.datasetFileOverride to the correct dataset file.']);
end

fprintf('Model was trained using dataset:\n%s\n', datasetFile);

%% ======================== LOAD NORMALIZATION ONCE =========================

[mu, sig, FsModel, Cmodel, Tmodel] = load_dataset_norm_and_size(datasetFile);

fprintf('\nModel/dataset settings:\n');
fprintf('FsModel = %.3f Hz\n', FsModel);
fprintf('Cmodel  = %d channels\n', Cmodel);
fprintf('Tmodel  = %d time samples\n', Tmodel);

nPre  = round(abs(opts.modelWinSec(1)) * FsModel);
nPost = Tmodel - nPre - 1;

if nPost < 0
    error('Bad modelWinSec/Tmodel combination. nPost < 0.');
end

fprintf('Model window: %.3f to %.3f sec | nPre=%d | nPost=%d | total=%d\n', ...
    opts.modelWinSec(1), opts.modelWinSec(2), nPre, nPost, nPre+nPost+1);

if (nPre + nPost + 1) ~= Tmodel
    warning('Window sample count does not exactly match Tmodel. Check opts.modelWinSec.');
end

mu = single(mu(:));
sig = single(sig(:));
sig(sig == 0 | isnan(sig)) = 1;
mu(isnan(mu)) = 0;

%% ======================== BATCH PROCESS =================================

Log = {};

for ff = 1:numel(rawFiles)

    rawFile = rawFiles{ff};
    [~, rawBase, ~] = fileparts(rawFile);

    fprintf('\n============================================================\n');
    fprintf('[%d/%d] Processing:\n%s\n', ff, numel(rawFiles), rawFile);
    fprintf('============================================================\n');

    matFile = fullfile(opts.outDir, [rawBase sprintf('_model_candidates_gt%d_refined_full.mat', round(opts.probThreshold*100))]);

    if opts.skipExisting && isfile(matFile)
        fprintf('Output exists, skipping:\n%s\n', matFile);
        Log(end+1,:) = {rawFile, 'SKIPPED_EXISTS', NaN, NaN, matFile}; %#ok<SAGROW>
        continue;
    end

    try
        [nEvents, nWindows, outFiles] = process_one_fif( ...
            rawFile, net, mu, sig, FsModel, Cmodel, Tmodel, nPre, nPost, opts);

        Log(end+1,:) = {rawFile, 'OK', nWindows, nEvents, outFiles.matFile}; %#ok<SAGROW>

    catch ME
        warning('Failed on file:\n%s\n%s', rawFile, ME.message);
        Log(end+1,:) = {rawFile, ['FAILED: ' ME.message], NaN, NaN, ''}; %#ok<SAGROW>
        continue;
    end
end

%% ======================== SAVE BATCH LOG =================================

Tlog = cell2table(Log, 'VariableNames', ...
    {'RawFile','Status','NWindows','NEvents','MatFile'});

logFile = fullfile(opts.outDir, ['Step6_batch_log_' datestr(now,'yyyymmdd_HHMMSS') '.csv']);
writetable(Tlog, logFile);

fprintf('\nBatch complete.\n');
fprintf('Saved log:\n%s\n', logFile);

disp(Tlog);

%% ========================================================================
% Main per-file processor
% ========================================================================

function [nEvents, nWin, outFiles] = process_one_fif( ...
    rawFile, net, mu, sig, FsModel, Cmodel, Tmodel, nPre, nPost, opts)

    [~, rawBase, ~] = fileparts(rawFile);

    %% ---- Load continuous MEG ----
    fprintf('\nLoading continuous MEG:\n%s\n', rawFile);

    [dataMat, chanLabels, FsData] = load_continuous_meg(rawFile, opts, FsModel);

    fprintf('Loaded MEG: %d channels x %d samples | Fs=%.3f Hz\n', ...
        size(dataMat,1), size(dataMat,2), FsData);

    if size(dataMat,1) < Cmodel
        error('Data has %d channels but model expects %d.', size(dataMat,1), Cmodel);
    end

    dataMat = dataMat(1:Cmodel, :);
    chanLabels = chanLabels(1:Cmodel);

    if opts.useTimeRange
        s1 = max(1, round(opts.timeRangeSec(1) * FsModel) + 1);
        s2 = min(size(dataMat,2), round(opts.timeRangeSec(2) * FsModel));
        dataMat = dataMat(:, s1:s2);
        timeOffsetSec = (s1 - 1) / FsModel;
    else
        timeOffsetSec = 0;
    end

    nData = size(dataMat,2);

    %% ---- Sliding-window prediction ----
    stepSamp = max(1, round(opts.stepSec * FsModel));

    centerSamples = (nPre+1):stepSamp:(nData-nPost);
    nWin = numel(centerSamples);

    fprintf('\nSliding prediction:\n');
    fprintf('Number of windows: %d\n', nWin);
    fprintf('Step: %.3f sec = %d samples\n', opts.stepSec, stepSamp);

    probSpike = zeros(nWin,1,'single');
    predLabel = zeros(nWin,1);

    for s = 1:opts.miniBatchWindows:nWin

        e = min(s + opts.miniBatchWindows - 1, nWin);
        idxWin = s:e;

        Xbatch = zeros(numel(idxWin), Cmodel, Tmodel, 'single');

        for ii = 1:numel(idxWin)

            c = centerSamples(idxWin(ii));
            idx = (c-nPre):(c+nPost);

            xi = single(dataMat(:, idx));
            xi = bsxfun(@rdivide, bsxfun(@minus, xi, mu), sig);

            Xbatch(ii,:,:) = reshape(xi, [1 Cmodel Tmodel]);
        end

        dlX = make_dlX_for_model(net, Xbatch);

        dlY = predict(net, dlX);
        P = gather(extractdata(dlY));

        probSpike(idxWin) = single(P(2,:)');

        [~, pc] = max(P, [], 1);
        predLabel(idxWin) = pc(:) - 1;

        if mod(e, 5000) == 0 || e == nWin
            fprintf('  predicted %d / %d windows\n', e, nWin);
        end
    end

    centerTimeSec = double(centerSamples(:)-1) ./ FsModel + timeOffsetSec;

    %% ---- Pick strong probability peaks ----
    fprintf('\nSelecting strong probability peaks...\n');

    smoothWinPts = max(1, round(opts.probSmoothSec / opts.stepSec));
    probForPick = movmean(double(probSpike), smoothWinPts);

    idxPeaksAll = find_local_maxima(probForPick);
    idxStrongPeaks = idxPeaksAll(probSpike(idxPeaksAll) >= opts.probThreshold);

    idxStrongAll = find(probSpike >= opts.probThreshold);

    fprintf('All windows >= %.2f: %d\n', opts.probThreshold, numel(idxStrongAll));
    fprintf('Strong local probability peaks: %d\n', numel(idxStrongPeaks));

    if isempty(idxStrongPeaks)
        warning('No strong local peaks found. Falling back to all above-threshold windows.');
        idxCand = idxStrongAll(:);
    else
        idxCand = idxStrongPeaks(:);
    end

    minSepWin = max(1, round(opts.minEventSepSec / opts.stepSec));
    idxCand = nonmax_suppress_idx(idxCand, double(probSpike), minSepWin);

    [~, ordTime] = sort(centerTimeSec(idxCand), 'ascend');
    idxCand = idxCand(ordTime);

    fprintf('Candidates after probability NMS: %d\n', numel(idxCand));

    %% ---- Refine to local GFP peak ----
    if opts.refineToGfpPeak && ~isempty(idxCand)

        fprintf('\nRefining candidate times to local GFP peak...\n');

        refineWinSamp = round(opts.refineWinSec * FsModel);

        refinedSamplesLocal = zeros(numel(idxCand),1);
        refinedTimes = zeros(numel(idxCand),1);

        for ii = 1:numel(idxCand)

            c = centerSamples(idxCand(ii));

            lo = max(1, c - refineWinSamp);
            hi = min(size(dataMat,2), c + refineWinSamp);

            seg = single(dataMat(:, lo:hi));
            seg = bsxfun(@rdivide, bsxfun(@minus, seg, mu), sig);

            gfp = sqrt(mean(seg.^2, 1));

            [~, imax] = max(gfp);

            refinedSamplesLocal(ii) = lo + imax - 1;
            refinedTimes(ii) = double(refinedSamplesLocal(ii)-1) ./ FsModel + timeOffsetSec;
        end

    else
        refinedSamplesLocal = centerSamples(idxCand(:))';
        refinedTimes = centerTimeSec(idxCand(:));
    end

    %% ---- Deduplicate after GFP refinement ----
    if ~isempty(idxCand)

        fprintf('\nRemoving duplicate refined events...\n');

        dedupMinSepSamp = round(opts.refineDedupSec * FsModel);
        probCand = double(probSpike(idxCand));

        [~, ordProb] = sort(probCand, 'descend');

        keep = false(numel(idxCand),1);
        keptSamples = [];

        for kk = 1:numel(ordProb)

            ii = ordProb(kk);
            s  = refinedSamplesLocal(ii);

            if isempty(keptSamples) || all(abs(s - keptSamples) > dedupMinSepSamp)
                keep(ii) = true;
                keptSamples(end+1,1) = s; %#ok<AGROW>
            end
        end

        idxCand = idxCand(keep);
        refinedSamplesLocal = refinedSamplesLocal(keep);
        refinedTimes = refinedTimes(keep);
        probCand = probCand(keep);

        [refinedTimes, ordTime] = sort(refinedTimes, 'ascend');
        idxCand = idxCand(ordTime);
        refinedSamplesLocal = refinedSamplesLocal(ordTime);
        probCand = probCand(ordTime);

    else
        probCand = [];
        refinedSamplesLocal = [];
        refinedTimes = [];
    end

    nEvents = numel(idxCand);

    fprintf('Final refined/deduplicated events exported: %d\n', nEvents);

    %% ---- Header/event coordinate conversion ----
    hdr = ft_read_header(rawFile);

    if isfield(hdr, 'orig') && isfield(hdr.orig, 'sfreq')
        fs_hdr = hdr.orig.sfreq;
    else
        fs_hdr = hdr.Fs;
    end

    try
        first_samp = double(hdr.orig.raw.first_samp);
    catch
        first_samp = 0;
    end

    candTimes = refinedTimes(:);
    candSamps = floor(candTimes .* fs_hdr) + first_samp;
    candTsecs = candTimes + first_samp ./ fs_hdr;

    %% ---- Save outputs ----
    tag = sprintf('gt%d_refined', round(opts.probThreshold*100));

    csvFile  = fullfile(opts.outDir, [rawBase '_model_candidates_' tag '_timepoints.csv']);
    xlsxFile = fullfile(opts.outDir, [rawBase '_model_candidates_' tag '_timepoints.xlsx']);
    eveFile  = fullfile(opts.outDir, [rawBase '_model_candidates_' tag '_spk_DS.eve']);
    matFile  = fullfile(opts.outDir, [rawBase '_model_candidates_' tag '_full.mat']);

    Treview = table();
    Treview.sample = round(candSamps(:));
    Treview.time_sec = candTsecs(:);

    writetable(Treview, csvFile);
    writetable(Treview, xlsxFile);

    fprintf('\nSaved minimal review CSV:\n%s\n', csvFile);
    fprintf('Saved minimal review Excel:\n%s\n', xlsxFile);

    % mbrowse-compatible event file
    fid = fopen(eveFile, 'w');
    assert(fid > 0, 'Cannot open %s for writing.', eveFile);

    fprintf(fid, '%d\t%f\t%d\t%d\t%s\n', ...
        first_samp, first_samp/fs_hdr, 0, 0, 'test');

    for ii = 1:numel(candSamps)
        fprintf(fid, '%d\t%f\t%d\t%d\n', ...
            round(candSamps(ii)), candTsecs(ii), 0, opts.eventCode);
    end

    fclose(fid);

    fprintf('Saved mbrowse-compatible event file:\n%s\n', eveFile);

    Tfull = table();
    Tfull.rawFile = repmat({rawFile}, numel(idxCand), 1);
    Tfull.rawBase = repmat({rawBase}, numel(idxCand), 1);
    Tfull.modelCenterSampleLocal = centerSamples(idxCand)';
    Tfull.modelCenterTimeSec = centerTimeSec(idxCand);
    Tfull.refinedSampleLocal = refinedSamplesLocal(:);
    Tfull.eventSample = round(candSamps(:));
    Tfull.eventTimeSec = candTsecs(:);
    Tfull.probSpike = double(probSpike(idxCand));
    Tfull.predictedLabel = predLabel(idxCand);
    Tfull.reviewLabel = nan(numel(idxCand),1);
    Tfull.reviewNotes = repmat({''}, numel(idxCand), 1);

    if opts.saveFullProbabilityTrace
        Tprob = table();
        Tprob.centerSampleLocal = centerSamples(:);
        Tprob.centerTimeSec = centerTimeSec(:);
        Tprob.probSpike = double(probSpike(:));
        Tprob.predictedLabel = predLabel(:);
    else
        Tprob = table();
    end

    save(matFile, ...
        'Treview','Tfull','Tprob','probSpike','predLabel','centerSamples','centerTimeSec', ...
        'idxCand','refinedSamplesLocal','refinedTimes','probCand', ...
        'eveFile','csvFile','xlsxFile', ...
        'opts','chanLabels','FsModel','Cmodel','Tmodel','nPre','nPost', ...
        '-v7.3');

    fprintf('Saved full MAT for feedback/retraining:\n%s\n', matFile);

    % Quick plot per file
    fig = figure('Visible','off');
    plot(centerTimeSec, probSpike);
    xlabel('Time (s)');
    ylabel('Spike probability');
    title(sprintf('Model spike probability: %s', rawBase), 'Interpreter','none');
    grid on;
    hold on;
    yline(opts.probThreshold, '--');
    if ~isempty(idxCand)
        scatter(refinedTimes, probSpike(idxCand), 25, 'filled');
    end
    hold off;

    pngFile = fullfile(opts.outDir, [rawBase '_model_candidates_' tag '_probability.png']);
    saveas(fig, pngFile);
    close(fig);

    outFiles = struct();
    outFiles.csvFile = csvFile;
    outFiles.xlsxFile = xlsxFile;
    outFiles.eveFile = eveFile;
    outFiles.matFile = matFile;
    outFiles.pngFile = pngFile;
end

%% ========================================================================
% Helper functions
% ========================================================================

function [mu, sig, Fs, C, T] = load_dataset_norm_and_size(datasetFile)

    info = whos('-file', datasetFile);
    names = {info.name};

    if all(ismember({'X_train','mu','sig','Fs'}, names))

        m = matfile(datasetFile);
        sz = size(m, 'X_train');

        C = sz(2);
        T = sz(3);

        S = load(datasetFile, 'mu','sig','Fs');
        mu = S.mu;
        sig = S.sig;
        Fs = S.Fs;

    elseif ismember('dataset', names)

        S = load(datasetFile, 'dataset');
        dataset = S.dataset;

        if ~isfield(dataset,'X_train')
            error('dataset does not contain X_train.');
        end

        sz = size(dataset.X_train);

        C = sz(2);
        T = sz(3);

        mu = dataset.mu;
        sig = dataset.sig;
        Fs = dataset.Fs;

    else
        error('Cannot determine dataset format for: %s', datasetFile);
    end
end


function [x, chanLabels, Fs] = load_continuous_meg(rawFile, opts, targetFs)

    cfg = [];
    cfg.dataset = rawFile;
    cfg.continuous = 'yes';
    cfg.channel = 'MEG';

    % Match training preprocessing if model was trained on bandpass data
    if isfield(opts, 'useBandpass') && opts.useBandpass
        fprintf('Applying bandpass filter %.1f-%.1f Hz before prediction...\n', ...
            opts.bpFreq(1), opts.bpFreq(2));

        cfg.bpfilter   = 'yes';
        cfg.bpfreq     = opts.bpFreq;
        cfg.bpfilttype = 'but';
        cfg.bpfiltord  = opts.bpOrder;
        cfg.bpfiltdir  = 'twopass';
    end

    data = ft_preprocessing(cfg);

    Fs = data.fsample;

    if abs(Fs - targetFs) > 1e-6
        fprintf('Resampling from %.3f Hz to %.3f Hz...\n', Fs, targetFs);

        cfg = [];
        cfg.resamplefs = targetFs;
        cfg.detrend = 'no';
        cfg.demean = 'no';

        data = ft_resampledata(cfg, data);
        Fs = data.fsample;
    end

    x = data.trial{1};
    chanLabels = data.label(:);

    if isempty(x)
        error('Loaded empty MEG matrix from %s', rawFile);
    end
end

% function [x, chanLabels, Fs] = load_continuous_meg(rawFile, opts, targetFs)
% 
%     cfg = [];
%     cfg.dataset = rawFile;
%     cfg.continuous = 'yes';
%     cfg.channel = 'MEG';
% 
%     data = ft_preprocessing(cfg);
% 
%     Fs = data.fsample;
% 
%     if abs(Fs - targetFs) > 1e-6
%         fprintf('Resampling from %.3f Hz to %.3f Hz...\n', Fs, targetFs);
% 
%         cfg = [];
%         cfg.resamplefs = targetFs;
%         cfg.detrend = 'no';
%         cfg.demean = 'no';
% 
%         data = ft_resampledata(cfg, data);
%         Fs = data.fsample;
%     end
% 
%     x = data.trial{1};
%     chanLabels = data.label(:);
% 
%     if isempty(x)
%         error('Loaded empty MEG matrix from %s', rawFile);
%     end
% end

function dlX = make_dlX_for_model(net, Xbatch)

    inputLayer = net.Layers(1);
    inputClass = class(inputLayer);

    B = size(Xbatch,1);
    C = size(Xbatch,2);
    T = size(Xbatch,3);

    persistent alreadyPrinted
    if isempty(alreadyPrinted)
        fprintf('Detected model input layer: %s\n', inputClass);
        alreadyPrinted = true;
    end

    if contains(inputClass, 'ImageInputLayer')

        X = reshape(Xbatch, [B C T 1]);
        X = permute(X, [2 3 4 1]);
        dlX = dlarray(single(X), 'SSCB');

    elseif contains(inputClass, 'SequenceInputLayer')

        X = permute(Xbatch, [2 3 1]);
        dlX = dlarray(single(X), 'CTB');

    else
        error('Unsupported input layer type: %s', inputClass);
    end
end

function idx = find_local_maxima(x)

    x = double(x(:));

    if numel(x) < 3
        idx = [];
        return;
    end

    idx = find(x(2:end-1) >= x(1:end-2) & x(2:end-1) >= x(3:end)) + 1;
end

function idxKeep = nonmax_suppress_idx(idxCand, scores, minSepIdx)

    idxCand = idxCand(:);

    if isempty(idxCand)
        idxKeep = idxCand;
        return;
    end

    scores = double(scores(:));

    [~, ord] = sort(scores(idxCand), 'descend');

    keep = false(numel(idxCand),1);
    kept = [];

    for kk = 1:numel(ord)

        ii = ord(kk);
        c = idxCand(ii);

        if isempty(kept) || all(abs(c - kept) > minSepIdx)
            keep(ii) = true;
            kept(end+1,1) = c; %#ok<AGROW>
        end
    end

    idxKeep = idxCand(keep);
end