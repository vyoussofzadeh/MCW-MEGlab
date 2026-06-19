%% ------------------------------------------------------------------------
% Step7b_Create_Feedback_Epochs_From_Reviewed_Events.m
%
% Convert Step7a interactive review output into model-shaped feedback epochs.
%
% Input:
%   Step7a output:
%       *_interactive_reviewed.mat
%
% Expected Step7a MAT contents:
%   Treviewed        table with reviewLabel
%   eventSampleLocal local sample index for each reviewed candidate
%   opts6            Step 6 options, usually including rawFile/modelFile
%   rawFile          raw FIF path, optional fallback
%
% Output:
%   feedback_round01_epochs.mat
%
% Output variables:
%   feedback.X       N x C x T single, normalized like the training set
%   feedback.y       N x 1 double, 1 spike / 0 no-spike
%   feedback.T       reviewed rows used for extraction
%   X_feedback       duplicate of feedback.X for older Step8 scripts
%   y_feedback       duplicate of feedback.y for older Step8 scripts
% -------------------------------------------------------------------------

clc; clear;

%% ======================== USER SETTINGS =================================

opts.ftRoot  = '/home/vyoussofzadeh/Desktop/research_workspace/tools/fieldtrip/fieldtrip_2022/';
opts.mneRoot = '/MEG_data/MEG_Tools/mne';

% Folder containing Step7a *_interactive_reviewed.mat files.
opts.reviewDir = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data/';
opts.reviewPattern = '*_interactive_reviewed.mat';

% Output feedback MAT file.
opts.outFile = fullfile(opts.reviewDir, 'feedback_round01_epochs.mat');

% Usually leave these empty. Use them when paths in the MAT files point to
% another machine/location.
opts.rawFileOverride     = '';
opts.modelFileOverride   = '';
opts.datasetFileOverride = '';

% Match Step 6 display/model preprocessing before extracting feedback epochs.
opts.matchStep6Preprocessing = true;

% Ignore unsure/unreviewed labels.
opts.keepLabels = [0 1];

% Save a small CSV summary beside opts.outFile.
opts.saveQuickSummaryCsv = true;

%% ======================== INITIALIZE ====================================

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

outDir = fileparts(opts.outFile);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ======================== FIND REVIEW FILES =============================

L = dir(fullfile(opts.reviewDir, opts.reviewPattern));

if isempty(L)
    error('No Step7a review files found:\n%s', fullfile(opts.reviewDir, opts.reviewPattern));
end

reviewFiles = fullfile({L.folder}, {L.name})';

fprintf('\nFound %d Step7a reviewed file(s):\n', numel(reviewFiles));
for ii = 1:numel(reviewFiles)
    fprintf('[%03d] %s\n', ii, reviewFiles{ii});
end

%% ======================== LOAD FIRST FILE METADATA =======================

Sfirst = load(reviewFiles{1});

if ~isfield(Sfirst, 'opts6')
    error('Step7a MAT file must contain opts6: %s', reviewFiles{1});
end

opts6 = Sfirst.opts6;

rawFile = determine_raw_file(Sfirst, opts);
modelFile = determine_model_file(opts6, opts);
datasetFile = determine_dataset_file(modelFile, opts);

fprintf('\nRaw file:\n%s\n', rawFile);
fprintf('\nModel file:\n%s\n', modelFile);
fprintf('\nDataset file:\n%s\n', datasetFile);

%% ======================== MODEL SIZE / NORMALIZATION =====================

[mu, sig, FsModel, Cmodel, Tmodel] = load_dataset_norm_and_size(datasetFile);

mu = single(mu(:));
sig = single(sig(:));
sig(sig == 0 | isnan(sig)) = 1;
mu(isnan(mu)) = 0;

[nPre, nPost] = determine_epoch_window(Sfirst, opts6, FsModel, Tmodel);

fprintf('\nFeedback epoch target:\n');
fprintf('FsModel = %.3f Hz\n', FsModel);
fprintf('Cmodel  = %d channels\n', Cmodel);
fprintf('Tmodel  = %d samples\n', Tmodel);
fprintf('nPre    = %d samples\n', nPre);
fprintf('nPost   = %d samples\n', nPost);

%% ======================== LOAD RAW MEG ONCE ==============================

fprintf('\nLoading raw MEG for feedback extraction...\n%s\n', rawFile);

[dataMat, chanLabels, FsData] = load_continuous_meg(rawFile, FsModel, opts6, opts);

if size(dataMat,1) < Cmodel
    error('Raw data has %d MEG channels but model expects %d.', size(dataMat,1), Cmodel);
end

dataMat = single(dataMat(1:Cmodel,:));
chanLabels = chanLabels(1:Cmodel);

fprintf('Loaded MEG: %d channels x %d samples | Fs=%.3f Hz\n', ...
    size(dataMat,1), size(dataMat,2), FsData);

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

%% ======================== EXTRACT ALL REVIEWED EPOCHS ===================

X_cells = {};
y_cells = {};
T_cells = {};
center_cells = {};
source_cells = {};
skipTables = {};

for ff = 1:numel(reviewFiles)

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Processing review file %d/%d:\n%s\n', ff, numel(reviewFiles), reviewFiles{ff});

    Srev = load(reviewFiles{ff});
    [Treview, reviewSource] = load_step7a_review_table(Srev, reviewFiles{ff});

    reviewLabel = parse_review_labels(Treview.reviewLabel);
    keep = ismember(reviewLabel, opts.keepLabels);

    fprintf('Rows total: %d | usable labels: %d | ignored: %d\n', ...
        height(Treview), sum(keep), sum(~keep));

    Treview = Treview(keep,:);
    reviewLabel = reviewLabel(keep);

    if height(Treview) == 0
        continue;
    end

    [centerSampleLocal, sampleSource] = get_center_samples_for_review( ...
        Srev, Treview, keep, FsModel, fs_hdr, first_samp);

    [X_feedback_i, y_feedback_i, Tfeedback_i, center_i, skip_i] = extract_feedback_epochs( ...
        Treview, reviewLabel, centerSampleLocal, dataMat, mu, sig, nPre, nPost, Cmodel, Tmodel);

    Tfeedback_i.reviewFile = repmat(string(reviewFiles{ff}), height(Tfeedback_i), 1);
    Tfeedback_i.reviewSource = repmat(string(reviewSource), height(Tfeedback_i), 1);
    Tfeedback_i.sampleMappingSource = repmat(string(sampleSource), height(Tfeedback_i), 1);

    fprintf('Extracted: %d | spike=%d | no-spike=%d | boundary skipped=%d | bad sample=%d\n', ...
        size(X_feedback_i,1), sum(y_feedback_i==1), sum(y_feedback_i==0), ...
        sum(skip_i.reason == "boundary"), sum(skip_i.reason == "bad_sample"));

    X_cells{end+1,1} = X_feedback_i; %#ok<SAGROW>
    y_cells{end+1,1} = y_feedback_i; %#ok<SAGROW>
    T_cells{end+1,1} = Tfeedback_i; %#ok<SAGROW>
    center_cells{end+1,1} = center_i; %#ok<SAGROW>
    source_cells{end+1,1} = repmat(string(sampleSource), numel(center_i), 1); %#ok<SAGROW>
    skip_i.reviewFile = repmat(string(reviewFiles{ff}), height(skip_i), 1);
    skipTables{end+1,1} = skip_i; %#ok<SAGROW>
end

if isempty(X_cells)
    error('No usable reviewed events were found.');
end

X_feedback = cat(1, X_cells{:});
y_feedback = cat(1, y_cells{:});
Tfeedback = vertcat(T_cells{:});
centerSampleLocal = vertcat(center_cells{:});
sampleMappingSource = vertcat(source_cells{:});

if isempty(skipTables)
    Tskip = table();
else
    Tskip = vertcat(skipTables{:});
end

if size(X_feedback,1) == 0
    error('No valid feedback epochs extracted.');
end

fprintf('\n============================================================\n');
fprintf('Feedback epochs extracted:\n');
fprintf('Total extracted: %d\n', size(X_feedback,1));
fprintf('True spikes:     %d\n', sum(y_feedback==1));
fprintf('No-spike:        %d\n', sum(y_feedback==0));

if numel(unique(y_feedback)) < 2
    warning('Feedback contains only one class. Fine-tuning can become biased.');
end

%% ======================== SAVE FEEDBACK =================================

feedback = struct();
feedback.X = X_feedback;
feedback.y = y_feedback;
feedback.T = Tfeedback;

feedback.centerSampleLocal = centerSampleLocal;
feedback.sampleMappingSource = sampleMappingSource;
feedback.skipped = Tskip;

feedback.rawFile = rawFile;
feedback.modelFile = modelFile;
feedback.datasetFile = datasetFile;
feedback.reviewFiles = reviewFiles;

feedback.Fs = FsModel;
feedback.C = Cmodel;
feedback.Tsamples = Tmodel;
feedback.nPre = nPre;
feedback.nPost = nPost;
feedback.chanLabels = chanLabels;

feedback.mu = mu;
feedback.sig = sig;
feedback.opts_step7b = opts;
feedback.opts_step6 = opts6;

feedback.matchStep6Preprocessing = opts.matchStep6Preprocessing;
feedback.useBandpass = isfield(opts6, 'useBandpass') && opts6.useBandpass;

if feedback.useBandpass
    feedback.bpFreq = opts6.bpFreq;
    if isfield(opts6, 'bpOrder')
        feedback.bpOrder = opts6.bpOrder;
    else
        feedback.bpOrder = NaN;
    end
else
    feedback.bpFreq = [];
    feedback.bpOrder = [];
end

feedback.created = datestr(now);

save(opts.outFile, ...
    'feedback','X_feedback','y_feedback','Tfeedback', ...
    'centerSampleLocal','sampleMappingSource','Tskip', ...
    'rawFile','modelFile','datasetFile', ...
    'FsModel','Cmodel','Tmodel','nPre','nPost','chanLabels', ...
    'opts','opts6', ...
    '-v7.3');

fprintf('\nSaved feedback epochs:\n%s\n', opts.outFile);

%% ---- Optional quick summary CSV ----------------------------------------

Tsummary = table();
Tsummary.N_feedback = size(X_feedback,1);
Tsummary.N_spike = sum(y_feedback==1);
Tsummary.N_nospike = sum(y_feedback==0);
Tsummary.Fs = FsModel;
Tsummary.C = Cmodel;
Tsummary.T = Tmodel;
Tsummary.nPre = nPre;
Tsummary.nPost = nPost;
Tsummary.N_skipped = height(Tskip);
Tsummary.UseBandpass = feedback.useBandpass;

if feedback.useBandpass
    Tsummary.BpLowHz = opts6.bpFreq(1);
    Tsummary.BpHighHz = opts6.bpFreq(2);
else
    Tsummary.BpLowHz = NaN;
    Tsummary.BpHighHz = NaN;
end

if opts.saveQuickSummaryCsv
    [p, b, ~] = fileparts(opts.outFile);
    summaryCsv = fullfile(p, [b '_summary.csv']);
    writetable(Tsummary, summaryCsv);
    fprintf('Saved feedback summary:\n%s\n', summaryCsv);
end

disp(Tsummary);

fprintf('\nStep 7b complete. Use opts.outFile as Step8 opts.feedbackFile.\n');

%% ========================================================================
% Helper functions
% ========================================================================

function rawFile = determine_raw_file(S, opts)

    rawFile = '';

    if ~isempty(opts.rawFileOverride)
        rawFile = opts.rawFileOverride;
    elseif isfield(S, 'rawFile') && ~isempty(S.rawFile)
        rawFile = S.rawFile;
    elseif isfield(S, 'opts6') && isfield(S.opts6, 'rawFile') && ~isempty(S.opts6.rawFile)
        rawFile = S.opts6.rawFile;
    elseif isfield(S, 'Treviewed') && ismember('rawFile', S.Treviewed.Properties.VariableNames) ...
            && height(S.Treviewed) > 0
        rawFile = S.Treviewed.rawFile{1};
    end

    rawFile = char(rawFile);

    if isempty(rawFile) || ~isfile(rawFile)
        error(['Could not determine/find raw FIF file.\n' ...
               'Set opts.rawFileOverride to the correct FIF path.']);
    end
end

function modelFile = determine_model_file(opts6, opts)

    modelFile = '';

    if ~isempty(opts.modelFileOverride)
        modelFile = opts.modelFileOverride;
    elseif isfield(opts6, 'modelFile') && ~isempty(opts6.modelFile)
        modelFile = opts6.modelFile;
    elseif isfield(opts6, 'trainedModelFile') && ~isempty(opts6.trainedModelFile)
        modelFile = opts6.trainedModelFile;
    end

    modelFile = char(modelFile);

    if isempty(modelFile) || ~isfile(modelFile)
        error(['Could not determine/find model file.\n' ...
               'Set opts.modelFileOverride to the correct Step 5 model path.']);
    end
end

function datasetFile = determine_dataset_file(modelFile, opts)

    if ~isempty(opts.datasetFileOverride)
        datasetFile = opts.datasetFileOverride;
    else
        M = load(modelFile, 'datasetFile');
        if isfield(M, 'datasetFile')
            datasetFile = M.datasetFile;
        else
            datasetFile = '';
        end
    end

    datasetFile = char(datasetFile);

    if isempty(datasetFile) || ~isfile(datasetFile)
        error(['Could not find datasetFile used for normalization.\n' ...
               'Set opts.datasetFileOverride to the correct dataset MAT file.']);
    end
end

function [T, source] = load_step7a_review_table(S, reviewFile)

    if isfield(S, 'Treviewed')
        T = S.Treviewed;
        source = 'MAT:Treviewed';
    elseif isfield(S, 'Tfeedback')
        T = S.Tfeedback;
        source = 'MAT:Tfeedback';
    elseif isfield(S, 'Tfull')
        T = S.Tfull;
        source = 'MAT:Tfull';
    else
        error('Review file does not contain Treviewed/Tfeedback/Tfull: %s', reviewFile);
    end

    if ~istable(T)
        error('Review object is not a table: %s', reviewFile);
    end

    if ~ismember('reviewLabel', T.Properties.VariableNames)
        error('Review table must contain reviewLabel: %s', reviewFile);
    end
end

function [nPre, nPost] = determine_epoch_window(S, opts6, FsModel, Tmodel)

    if isfield(S, 'nPre') && isfield(S, 'nPost')
        nPre = round(double(S.nPre));
        nPost = round(double(S.nPost));
    elseif isfield(opts6, 'nPre') && isfield(opts6, 'nPost')
        nPre = round(double(opts6.nPre));
        nPost = round(double(opts6.nPost));
    elseif isfield(opts6, 'modelWinSec')
        modelWinSec = double(opts6.modelWinSec);
        nPre = round(abs(modelWinSec(1)) * FsModel);
        nPost = Tmodel - nPre - 1;
    else
        error('Could not determine nPre/nPost from Step7a/Step6 metadata.');
    end

    if (nPre + nPost + 1) ~= Tmodel
        error('Bad feedback window: nPre+nPost+1 = %d, but model expects Tmodel = %d.', ...
            nPre+nPost+1, Tmodel);
    end
end

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

        if ~isfield(dataset, 'X_train')
            error('dataset exists but does not contain X_train.');
        end

        sz = size(dataset.X_train);

        C = sz(2);
        T = sz(3);

        mu = dataset.mu;
        sig = dataset.sig;
        Fs = dataset.Fs;

    else
        error('Cannot determine dataset format: %s', datasetFile);
    end
end

function [x, chanLabels, Fs] = load_continuous_meg(rawFile, targetFs, opts6, opts7)

    cfg = [];
    cfg.dataset = rawFile;
    cfg.continuous = 'yes';
    cfg.channel = 'MEG';

    if opts7.matchStep6Preprocessing && isfield(opts6, 'useBandpass') && opts6.useBandpass

        if ~isfield(opts6, 'bpFreq')
            error('opts6.useBandpass=true but opts6.bpFreq is missing.');
        end

        fprintf('Applying feedback extraction bandpass %.1f-%.1f Hz...\n', ...
            opts6.bpFreq(1), opts6.bpFreq(2));

        cfg.bpfilter   = 'yes';
        cfg.bpfreq     = opts6.bpFreq;
        cfg.bpfilttype = 'but';

        if isfield(opts6, 'bpOrder')
            cfg.bpfiltord = opts6.bpOrder;
        else
            cfg.bpfiltord = 4;
        end

        cfg.bpfiltdir = 'twopass';
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
        error('Loaded empty MEG matrix from raw file.');
    end
end

function y = parse_review_labels(v)

    if isnumeric(v)
        y = double(v(:));
        return;
    end

    s = lower(strtrim(string(v)));
    y = nan(numel(s),1);

    pos = ismember(s, [ ...
        "1","yes","y","true","t", ...
        "spike","real","true spike","true_spike","tp"]);

    neg = ismember(s, [ ...
        "0","no","n","false","f", ...
        "nospike","no spike","no_spike", ...
        "artifact","false positive","false_positive","fp"]);

    y(pos) = 1;
    y(neg) = 0;

    nums = str2double(s);
    useNum = isnan(y) & ~isnan(nums);
    y(useNum) = nums(useNum);
end

function [centerSampleLocal, sourceName] = get_center_samples_for_review( ...
    S, Treview, keepMaskOriginal, Fs, fs_hdr, first_samp)

    if isfield(S, 'eventSampleLocal')
        allSamples = double(S.eventSampleLocal(:));
        if numel(allSamples) == numel(keepMaskOriginal)
            centerSampleLocal = round(allSamples(keepMaskOriginal));
            sourceName = 'Step7a:eventSampleLocal';
            return;
        end
    end

    [centerSampleLocal, sourceName] = get_candidate_local_samples_from_table( ...
        Treview, Fs, fs_hdr, first_samp);
end

function [X_feedback, y_feedback, Tfeedback, centerOut, Tskip] = extract_feedback_epochs( ...
    Treview, reviewLabel, centerSampleLocal, dataMat, mu, sig, nPre, nPost, Cmodel, Tmodel)

    Nreview = height(Treview);

    X_all = zeros(Nreview, Cmodel, Tmodel, 'single');
    valid = false(Nreview,1);
    skipReason = strings(Nreview,1);

    for ii = 1:Nreview

        c = centerSampleLocal(ii);

        if isnan(c) || ~isfinite(c)
            skipReason(ii) = "bad_sample";
            continue;
        end

        c = round(c);
        idx = (c-nPre):(c+nPost);

        if min(idx) < 1 || max(idx) > size(dataMat,2)
            skipReason(ii) = "boundary";
            continue;
        end

        xi = single(dataMat(:,idx));
        xi = bsxfun(@rdivide, bsxfun(@minus, xi, mu), sig);

        X_all(ii,:,:) = reshape(xi, [1 Cmodel Tmodel]);
        valid(ii) = true;
        skipReason(ii) = "";
    end

    X_feedback = X_all(valid,:,:);
    y_feedback = double(reviewLabel(valid));

    Tfeedback = Treview(valid,:);
    Tfeedback.feedbackLabel = y_feedback;
    Tfeedback.centerSampleLocal_feedback = centerSampleLocal(valid);

    centerOut = centerSampleLocal(valid);

    Tskip = table();
    Tskip.reviewRow = find(~valid);
    Tskip.reason = skipReason(~valid);
    Tskip.centerSampleLocal = centerSampleLocal(~valid);
end

function [centerSampleLocal, sourceName] = get_candidate_local_samples_from_table(T, Fs, fs_hdr, first_samp)

    name = find_var_name(T, { ...
        'refinedSampleLocal', ...
        'centerSampleLocal', ...
        'modelCenterSampleLocal'});

    if ~isempty(name)
        centerSampleLocal = round(col_to_double(T.(name)));
        sourceName = name;
        return;
    end

    name = find_var_name(T, { ...
        'eventSample', ...
        'sample', ...
        'samples', ...
        'event_sample'});

    if ~isempty(name)
        eventSampleAbs = col_to_double(T.(name));
        centerSampleLocal = round((eventSampleAbs - first_samp) ./ fs_hdr .* Fs) + 1;
        sourceName = name;
        return;
    end

    name = find_var_name(T, { ...
        'refinedTimeSec', ...
        'refinedTimes', ...
        'modelCenterTimeSec', ...
        'centerTimeSec'});

    if ~isempty(name)
        localTimeSec = col_to_double(T.(name));
        centerSampleLocal = round(localTimeSec .* Fs) + 1;
        sourceName = name;
        return;
    end

    name = find_var_name(T, { ...
        'eventTimeSec', ...
        'time_sec', ...
        'timesec', ...
        'time'});

    if ~isempty(name)
        eventTimeAbs = col_to_double(T.(name));
        eventTimeRel = eventTimeAbs - first_samp ./ fs_hdr;
        centerSampleLocal = round(eventTimeRel .* Fs) + 1;
        sourceName = name;
        return;
    end

    error('Could not find sample/time columns in review table.');
end

function name = find_var_name(T, aliases)

    vars = T.Properties.VariableNames;
    normVars = normalize_names(vars);
    normAliases = normalize_names(aliases);

    name = '';

    for ii = 1:numel(normAliases)
        idx = find(normVars == normAliases(ii), 1, 'first');
        if ~isempty(idx)
            name = vars{idx};
            return;
        end
    end
end

function out = normalize_names(names)

    names = string(names);
    out = lower(regexprep(names, '[^a-zA-Z0-9]', ''));
end

function x = col_to_double(v)

    if isnumeric(v)
        x = double(v(:));
    elseif iscell(v)
        x = str2double(string(v(:)));
    elseif isstring(v) || ischar(v) || iscategorical(v)
        x = str2double(string(v(:)));
    else
        x = double(v(:));
    end
end