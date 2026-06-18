%% ------------------------------------------------------------------------
% Step6_Run_Trained_Model_On_Continuous_MEG_Export_ReviewSheet.m
%
% Apply trained spike/no-spike model to one continuous MEG file.
%
% Outputs:
%   1) *_model_candidates_gt70_refined_timepoints.csv
%   2) *_model_candidates_gt70_refined_timepoints.xlsx
%   3) *_model_candidates_gt70_refined_spk_DS.eve
%   4) *_model_candidates_gt70_full.mat
%
% Excel/CSV contains only:
%   sample, time_sec
%
% .eve file follows the old mbrowse-compatible format:
%   first_samp    first_samp/fs_hdr    0    0    test
%   event_sample  event_time_sec       0    5555
% -------------------------------------------------------------------------

clc; clear;

%% ======================== USER SETTINGS =================================

opts.ftRoot  = '/home/vyoussofzadeh/Desktop/research_workspace/tools/fieldtrip/fieldtrip_2022/';
opts.mneRoot = '/MEG_data/MEG_Tools/mne';

% Trained model from Step 5
% opts.modelFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Squiggles/Step5_lowmem_2DCNN_20260518_170934.mat';
% opts.modelFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Squiggles/Step5_lowmem_CNN1D_20260527_153010.mat';
% opts.modelFile = '/home/vyoussofzadeh/Data/DL_model/Step5_lowmem_CNN1D_20260602_181546.mat';
% opts.modelFile = '/home/vyoussofzadeh/Data/DL_model/Step5_lowmem_CNN1D_20260603_210355.mat';
opts.modelFile = '/home/vyoussofzadeh/Data/DL_model/Step5_lowmem_CNN1D_bp5_50_m250_p500_20260604_134703.mat';


% If the model file points to a dataset path that does not exist on this machine,
% set this to the correct dataset used for training. Otherwise leave empty.
opts.datasetFileOverride = '';

% One continuous/sample FIF file to test
opts.rawFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data/baer_jessica_Run02_spont_eyesclosed_raw_t_sss_ecgClean_raw_DS.fif';

% Output folder
opts.outDir = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data';

% Sliding-window settings
opts.stepSec = 0.010;          % 50 ms step. Use 0.010 for finer timing.
% opts.modelWinSec = [-0.500 1.025];  % match the dataset used for model training
opts.modelWinSec = [-0.250 0.500];

% Strong candidate selection
opts.probThreshold = 0.70;     % only export candidates >= 70%
opts.probSmoothSec = 0.050;    % mild smoothing before finding probability peaks
opts.minEventSepSec = 0.300;   % keep only one probability peak within 300 ms

% Local peak refinement
opts.refineToGfpPeak = true;
opts.refineWinSec = 0.150;     % search +/-150 ms around model peak
opts.refineDedupSec = 0.100;   % after refinement, remove duplicates within 200 ms

% Event code for mbrowse/mne_browse_raw
opts.eventCode = 5555;

% Optional: process only part of recording for first test
opts.useTimeRange = false;
opts.timeRangeSec = [0 300];   % first 5 min if useTimeRange=true

% Save complete probability trace in MAT file
opts.saveFullProbabilityTrace = true;

%%
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

%% ======================== LOAD MODEL ====================================

fprintf('\nLoading model:\n%s\n', opts.modelFile);

M = load(opts.modelFile);

if isfield(M, 'bestNet')
    net = M.bestNet;
elseif isfield(M, 'net')
    net = M.net;
else
    error('Model file does not contain bestNet or net.');
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

%% ======================== LOAD NORMALIZATION =============================

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

%% ======================== LOAD CONTINUOUS MEG ============================

fprintf('\nLoading continuous MEG:\n%s\n', opts.rawFile);

[dataMat, chanLabels, FsData] = load_continuous_meg(opts.rawFile, opts, FsModel);

fprintf('Loaded MEG: %d channels x %d samples | Fs=%.3f Hz\n', ...
    size(dataMat,1), size(dataMat,2), FsData);

if size(dataMat,1) < Cmodel
    error('Data has %d channels but model expects %d.', size(dataMat,1), Cmodel);
end

% Keep first Cmodel MEG channels, matching training convention
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

%% ======================== SLIDING-WINDOW PREDICTION ======================

stepSamp = max(1, round(opts.stepSec * FsModel));

centerSamples = (nPre+1):stepSamp:(nData-nPost);
nWin = numel(centerSamples);

fprintf('\nSliding prediction:\n');
fprintf('Number of windows: %d\n', nWin);
fprintf('Step: %.3f sec = %d samples\n', opts.stepSec, stepSamp);

probSpike = zeros(nWin,1,'single');
predLabel = zeros(nWin,1);

miniBatchWindows = 128;

for s = 1:miniBatchWindows:nWin

    e = min(s + miniBatchWindows - 1, nWin);
    idxWin = s:e;

    Xbatch = zeros(numel(idxWin), Cmodel, Tmodel, 'single');

    for ii = 1:numel(idxWin)

        c = centerSamples(idxWin(ii));
        idx = (c-nPre):(c+nPost);

        xi = single(dataMat(:, idx));       % C x T
        xi = bsxfun(@rdivide, bsxfun(@minus, xi, mu), sig);

        Xbatch(ii,:,:) = reshape(xi, [1 Cmodel Tmodel]);
    end

    dlX = make_dlX_for_model(net, Xbatch);

    dlY = predict(net, dlX);
    P = gather(extractdata(dlY));   % 2 x B

    probSpike(idxWin) = single(P(2,:)');

    [~, pc] = max(P, [], 1);
    predLabel(idxWin) = pc(:) - 1;

    if mod(e, 5000) == 0 || e == nWin
        fprintf('  predicted %d / %d windows\n', e, nWin);
    end
end

centerTimeSec = double(centerSamples(:)-1) ./ FsModel + timeOffsetSec;

%% ======================== PICK STRONG PROBABILITY PEAKS ==================

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

% Keep only the strongest point within nearby groups
minSepWin = max(1, round(opts.minEventSepSec / opts.stepSec));
idxCand = nonmax_suppress_idx(idxCand, double(probSpike), minSepWin);

% Sort chronologically
[~, ordTime] = sort(centerTimeSec(idxCand), 'ascend');
idxCand = idxCand(ordTime);

fprintf('Candidates after probability NMS: %d\n', numel(idxCand));

%% ======================== REFINE CANDIDATES TO LOCAL GFP PEAK ============

if opts.refineToGfpPeak && ~isempty(idxCand)

    fprintf('\nRefining candidate times to local GFP peak...\n');

    refineWinSamp = round(opts.refineWinSec * FsModel);

    refinedSamplesLocal = zeros(numel(idxCand),1);
    refinedTimes = zeros(numel(idxCand),1);

    for ii = 1:numel(idxCand)

        c = centerSamples(idxCand(ii));  % local sample within dataMat

        lo = max(1, c - refineWinSamp);
        hi = min(size(dataMat,2), c + refineWinSamp);

        seg = single(dataMat(:, lo:hi));  % C x time

        % Same normalization as model
        seg = bsxfun(@rdivide, bsxfun(@minus, seg, mu), sig);

        % Global field power / global signal energy
        gfp = sqrt(mean(seg.^2, 1));

        [~, imax] = max(gfp);

        refinedSamplesLocal(ii) = lo + imax - 1;
        refinedTimes(ii) = double(refinedSamplesLocal(ii)-1) ./ FsModel + timeOffsetSec;
    end

else

    refinedSamplesLocal = centerSamples(idxCand(:))';
    refinedTimes = centerTimeSec(idxCand(:));
end

%% ======================== DEDUPLICATE AFTER GFP REFINEMENT ===============

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
            keptSamples(end+1,1) = s; 
        end
    end

    idxCand = idxCand(keep);
    refinedSamplesLocal = refinedSamplesLocal(keep);
    refinedTimes = refinedTimes(keep);
    probCand = probCand(keep);

    % Sort final events chronologically
    [refinedTimes, ordTime] = sort(refinedTimes, 'ascend');
    idxCand = idxCand(ordTime);
    refinedSamplesLocal = refinedSamplesLocal(ordTime);
    probCand = probCand(ordTime);

else
    probCand = [];
end

fprintf('Final refined/deduplicated events exported: %d\n', numel(idxCand));

%% ======================== HEADER INFO FOR EVENT EXPORT ===================

[~, rawBase, ~] = fileparts(opts.rawFile);

hdr = ft_read_header(opts.rawFile);

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

% Convert final refined event times to FIFF sample/time convention
candTimes = refinedTimes(:);
candSamps = floor(candTimes .* fs_hdr) + first_samp;
candTsecs = candTimes + first_samp ./ fs_hdr;

%% ======================== SAVE MINIMAL EXCEL/CSV =========================

% Excel/CSV only contains time point / sample
Treview = table();
Treview.sample = round(candSamps(:));
Treview.time_sec = candTsecs(:);

csvFile  = fullfile(opts.outDir, [rawBase '_model_candidates_gt70_refined_timepoints.csv']);
xlsxFile = fullfile(opts.outDir, [rawBase '_model_candidates_gt70_refined_timepoints.xlsx']);

writetable(Treview, csvFile);
writetable(Treview, xlsxFile);

fprintf('\nSaved minimal review CSV:\n%s\n', csvFile);
fprintf('Saved minimal review Excel:\n%s\n', xlsxFile);

%% ======================== SAVE MBROWSE/MNE EVENT FILE ====================

eveFile = fullfile(opts.outDir, [rawBase '_model_candidates_gt70_refined_spk_DS.eve']);

fid = fopen(eveFile, 'w');
assert(fid > 0, 'Cannot open %s for writing.', eveFile);

% Seed/header line for mbrowse compatibility
fprintf(fid, '%d\t%f\t%d\t%d\t%s\n', ...
    first_samp, first_samp/fs_hdr, 0, 0, 'test');

% Strong, refined, deduplicated candidate events only
for ii = 1:numel(candSamps)
    fprintf(fid, '%d\t%f\t%d\t%d\n', ...
        round(candSamps(ii)), candTsecs(ii), 0, opts.eventCode);
end

fclose(fid);

fprintf('Saved mbrowse-compatible refined event file:\n%s\n', eveFile);

%% ======================== SAVE FULL MAT FOR FEEDBACK =====================

matFile = fullfile(opts.outDir, [rawBase '_model_candidates_gt70_refined_full.mat']);

Tfull = table();
Tfull.rawFile = repmat({opts.rawFile}, numel(idxCand), 1);
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

%% ======================== QUICK PLOT =====================================

figure;
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

fprintf('\nStep 6 complete.\n');

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

function dlX = make_dlX_for_model(net, Xbatch)

    % Xbatch: B x C x T
    %
    % Supports:
    %   imageInputLayer    -> C x T x 1 x B, 'SSCB'
    %   sequenceInputLayer -> C x T x B,     'CTB'

    inputLayer = net.Layers(1);
    inputClass = class(inputLayer);

    B = size(Xbatch,1);
    C = size(Xbatch,2);
    T = size(Xbatch,3);

    fprintf_once_input_type(inputClass);

    if contains(inputClass, 'ImageInputLayer')

        % 2D CNN format
        X = reshape(Xbatch, [B C T 1]);   % B x C x T x 1
        X = permute(X, [2 3 4 1]);        % C x T x 1 x B
        dlX = dlarray(single(X), 'SSCB');

    elseif contains(inputClass, 'SequenceInputLayer')

        % CNN1D / LSTM format
        X = permute(Xbatch, [2 3 1]);     % C x T x B
        dlX = dlarray(single(X), 'CTB');

    else
        error('Unsupported input layer type: %s', inputClass);
    end
end

function fprintf_once_input_type(inputClass)
    persistent alreadyPrinted
    if isempty(alreadyPrinted)
        fprintf('Detected model input layer: %s\n', inputClass);
        alreadyPrinted = true;
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