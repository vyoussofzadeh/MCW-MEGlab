%% ------------------------------------------------------------------------
% Step7_Create_Feedback_Epochs_From_Reviewed_Events.m
%
% Read reviewed Step 6 candidate events, extract matching MEG epochs from
% the raw FIF file, normalize with the original training mu/sig, and save
% feedback data for Step 8 model fine-tuning/retraining.
%
% Input review file can be:
%   1) Excel/CSV from Step 6 with added reviewLabel column
%   2) Interactive review MAT file containing Treviewed
%
% Required review columns:
%   reviewLabel:
%       1 = true spike
%       0 = false positive / no-spike
%
% Event coordinate columns accepted:
%   refinedSampleLocal, modelCenterSampleLocal, centerSampleLocal
%   OR sample/eventSample
%   OR time_sec/eventTimeSec
%
% Output:
%   feedback_round01_epochs.mat
%
% Output variables:
%   feedback.X : N x channels x time
%   feedback.y : N x 1 labels, 1 spike / 0 no-spike
%   feedback.T : reviewed table used for extraction
% -------------------------------------------------------------------------

clc; clear;

%% ======================== USER SETTINGS =================================

opts.ftRoot  = '/home/vyoussofzadeh/Desktop/research_workspace/tools/fieldtrip/fieldtrip_2022/';
opts.mneRoot = '/MEG_data/MEG_Tools/mne';

% Reviewed file.
% Option A: manually edited Excel/CSV from Step 6 with reviewLabel column.
opts.reviewFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data/baer_jessica_Run02_spont_eyesclosed_raw_t_sss_ecgClean_raw_DS_model_candidates_gt70_refined_timepoints.xlsx';

% Option B: interactive review MAT file.
% opts.reviewFile = '/path/to/file_interactive_reviewed.mat';

% Full Step 6 MAT file. This contains opts, Tfull, nPre/nPost, etc.
opts.step6MatFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data/baer_jessica_Run02_spont_eyesclosed_raw_t_sss_ecgClean_raw_DS_model_candidates_gt70_refined_full.mat';

% Output feedback file
opts.outFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data/feedback_round01_epochs.mat';

% Overrides, usually leave empty.
% Use these if Step 6/model contains paths from another machine.
opts.rawFileOverride     = '';
opts.modelFileOverride   = '';
opts.datasetFileOverride = '';

% Match Step 6/model preprocessing.
% For BP model, this applies the same 5-50 Hz filtering before extracting
% feedback epochs.
opts.matchStep6Preprocessing = true;

% Save a small CSV summary next to the feedback MAT file.
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

if isempty(opts.outFile)
    [p, b, ~] = fileparts(opts.reviewFile);
    opts.outFile = fullfile(p, [b '_feedback_epochs.mat']);
end

outDir = fileparts(opts.outFile);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('\nChecking readers:\n');
which ft_preprocessing -all
which fiff_find_evoked -all

%% ======================== LOAD STEP 6 INFO ===============================

fprintf('\nLoading Step 6 MAT file:\n%s\n', opts.step6MatFile);

S6 = load(opts.step6MatFile);

if ~isfield(S6, 'opts')
    error('Step 6 MAT file does not contain opts.');
end

opts6 = S6.opts;

%% ---- Determine raw FIF file --------------------------------------------

rawFile = '';

if ~isempty(opts.rawFileOverride)
    rawFile = opts.rawFileOverride;

elseif isfield(opts6, 'rawFile') && ~isempty(opts6.rawFile)
    rawFile = opts6.rawFile;

elseif isfield(S6, 'Tfull') && ismember('rawFile', S6.Tfull.Properties.VariableNames) ...
        && height(S6.Tfull) > 0
    rawFile = S6.Tfull.rawFile{1};
end

if isempty(rawFile) || ~isfile(rawFile)
    error(['Could not determine/find raw FIF file.\n' ...
           'Set opts.rawFileOverride to the correct FIF path.']);
end

%% ---- Determine model file ----------------------------------------------

modelFile = '';

if ~isempty(opts.modelFileOverride)
    modelFile = opts.modelFileOverride;

elseif isfield(opts6, 'modelFile') && ~isempty(opts6.modelFile)
    modelFile = opts6.modelFile;
end

if isempty(modelFile) || ~isfile(modelFile)
    error(['Could not determine/find model file.\n' ...
           'Set opts.modelFileOverride to the correct Step 5 model path.']);
end

%% ---- Determine training dataset file -----------------------------------

M = load(modelFile, 'datasetFile');

if isfield(M, 'datasetFile')
    datasetFile = char(M.datasetFile);
else
    datasetFile = '';
end

if ~isempty(opts.datasetFileOverride)
    datasetFile = opts.datasetFileOverride;
end

if isempty(datasetFile) || ~isfile(datasetFile)
    error(['Could not find datasetFile used for normalization.\n' ...
           'Set opts.datasetFileOverride to the correct dataset MAT file.']);
end

fprintf('\nRaw file:\n%s\n', rawFile);
fprintf('\nModel file:\n%s\n', modelFile);
fprintf('\nDataset file:\n%s\n', datasetFile);

%% ======================== LOAD NORMALIZATION / SIZE ======================

[mu, sig, FsModel, Cmodel, Tmodel] = load_dataset_norm_and_size(datasetFile);

mu = single(mu(:));
sig = single(sig(:));
sig(sig == 0 | isnan(sig)) = 1;
mu(isnan(mu)) = 0;

%% ---- Get model window from Step 6 ---------------------------------------

if isfield(S6, 'nPre') && isfield(S6, 'nPost')
    nPre  = round(double(S6.nPre));
    nPost = round(double(S6.nPost));

elseif isfield(opts6, 'modelWinSec')
    modelWinSec = opts6.modelWinSec;
    nPre  = round(abs(modelWinSec(1)) * FsModel);
    nPost = Tmodel - nPre - 1;

else
    error('Could not determine nPre/nPost from Step 6 MAT file.');
end

if (nPre + nPost + 1) ~= Tmodel
    error('Bad feedback window: nPre+nPost+1 = %d, but Tmodel = %d.', ...
        nPre+nPost+1, Tmodel);
end

fprintf('\nFeedback epoch size:\n');
fprintf('FsModel = %.3f Hz\n', FsModel);
fprintf('Cmodel  = %d channels\n', Cmodel);
fprintf('Tmodel  = %d time samples\n', Tmodel);
fprintf('nPre    = %d samples\n', nPre);
fprintf('nPost   = %d samples\n', nPost);

if isfield(opts6, 'modelWinSec')
    expectedT = round((opts6.modelWinSec(2) - opts6.modelWinSec(1)) * FsModel) + 1;
    if expectedT ~= Tmodel
        warning('opts6.modelWinSec gives %d samples, but model expects Tmodel=%d.', ...
            expectedT, Tmodel);
    end
end

%% ======================== LOAD REVIEW TABLE ==============================

fprintf('\nLoading reviewed events:\n%s\n', opts.reviewFile);

[Treview, reviewSource] = load_review_table(opts.reviewFile);

fprintf('Review source: %s\n', reviewSource);
fprintf('Rows in review table: %d\n', height(Treview));

if isfield(S6, 'Tfull') && height(S6.Tfull) > 0 && height(Treview) ~= height(S6.Tfull)
    warning('Review table has %d rows but Step 6 Tfull has %d rows. Make sure these files match.', ...
        height(Treview), height(S6.Tfull));
end

labelName = find_var_name(Treview, {'reviewLabel','label','humanLabel','spike','feedbackLabel'});

if isempty(labelName)
    error('Review file must contain reviewLabel column: 1=true spike, 0=false positive/no-spike.');
end

reviewLabel = parse_review_labels(Treview.(labelName));

Treview.feedbackLabel = reviewLabel;

keepReview = reviewLabel == 0 | reviewLabel == 1;

Treview = Treview(keepReview,:);
reviewLabel = reviewLabel(keepReview);

fprintf('\nUsable reviewed events:\n');
fprintf('Total usable:      %d\n', height(Treview));
fprintf('True spikes:       %d\n', sum(reviewLabel==1));
fprintf('False positives:   %d\n', sum(reviewLabel==0));
fprintf('Ignored/unsure:    %d\n', sum(~keepReview));

if height(Treview) == 0
    error('No usable reviewed events found.');
end

%% ======================== LOAD RAW MEG ==================================

fprintf('\nLoading raw MEG for feedback extraction...\n%s\n', rawFile);

[dataMat, chanLabels, FsData] = load_continuous_meg(rawFile, FsModel, opts6, opts);

if size(dataMat,1) < Cmodel
    error('Raw data has %d channels but model expects %d.', size(dataMat,1), Cmodel);
end

dataMat = dataMat(1:Cmodel,:);
chanLabels = chanLabels(1:Cmodel);

fprintf('Loaded MEG: %d channels x %d samples | Fs=%.3f Hz\n', ...
    size(dataMat,1), size(dataMat,2), FsData);

%% ======================== MAP EVENTS TO LOCAL SAMPLES ====================

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

[centerSampleLocal, sampleMappingSource] = get_candidate_local_samples( ...
    Treview, FsModel, fs_hdr, first_samp);

fprintf('\nEvent sample mapping source: %s\n', sampleMappingSource);

%% ======================== EXTRACT FEEDBACK EPOCHS ========================

Nreview = height(Treview);

X_feedback_all = zeros(Nreview, Cmodel, Tmodel, 'single');
valid = false(Nreview,1);
skipReason = strings(Nreview,1);

for i = 1:Nreview

    c = centerSampleLocal(i);

    if isnan(c) || ~isfinite(c)
        skipReason(i) = "bad_sample";
        continue;
    end

    c = round(c);
    idx = (c-nPre):(c+nPost);

    if min(idx) < 1 || max(idx) > size(dataMat,2)
        skipReason(i) = "boundary";
        continue;
    end

    xi = single(dataMat(:,idx));  % C x T

    % Same normalization as model training.
    xi = bsxfun(@rdivide, bsxfun(@minus, xi, mu), sig);

    X_feedback_all(i,:,:) = reshape(xi, [1 Cmodel Tmodel]);
    valid(i) = true;
    skipReason(i) = "";
end

X_feedback = X_feedback_all(valid,:,:);
y_feedback = double(reviewLabel(valid));

Tfeedback = Treview(valid,:);
Tfeedback.feedbackLabel = y_feedback;
Tfeedback.centerSampleLocal_feedback = centerSampleLocal(valid);

centerSampleLocal = centerSampleLocal(valid);

fprintf('\nFeedback epochs extracted:\n');
fprintf('Total extracted:   %d\n', size(X_feedback,1));
fprintf('True spikes:       %d\n', sum(y_feedback==1));
fprintf('False positives:   %d\n', sum(y_feedback==0));
fprintf('Skipped boundary:  %d\n', sum(skipReason=="boundary"));
fprintf('Skipped bad sample:%d\n', sum(skipReason=="bad_sample"));

if size(X_feedback,1) == 0
    error('No valid feedback epochs extracted.');
end

%% ======================== SAVE FEEDBACK =================================

feedback = struct();

feedback.X = X_feedback;
feedback.y = y_feedback;
feedback.T = Tfeedback;

feedback.centerSampleLocal = centerSampleLocal;
feedback.sampleMappingSource = sampleMappingSource;

feedback.rawFile = rawFile;
feedback.modelFile = modelFile;
feedback.datasetFile = datasetFile;

feedback.Fs = FsModel;
feedback.C = Cmodel;
feedback.Tsamples = Tmodel;
feedback.nPre = nPre;
feedback.nPost = nPost;
feedback.chanLabels = chanLabels;

feedback.reviewFile = opts.reviewFile;
feedback.step6MatFile = opts.step6MatFile;
feedback.reviewSource = reviewSource;

feedback.opts_step7 = opts;
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
    'centerSampleLocal','rawFile','modelFile','datasetFile', ...
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

fprintf('\nStep 7 complete.\n');

%% ========================================================================
% Helper functions
% ========================================================================

function [T, source] = load_review_table(reviewFile)

    if ~isfile(reviewFile)
        error('Review file not found: %s', reviewFile);
    end

    [~, ~, ext] = fileparts(reviewFile);
    ext = lower(ext);

    switch ext

        case '.mat'

            S = load(reviewFile);

            if isfield(S, 'Treviewed')
                T = S.Treviewed;
                source = 'MAT:Treviewed';

            elseif isfield(S, 'Tfeedback')
                T = S.Tfeedback;
                source = 'MAT:Tfeedback';

            elseif isfield(S, 'Tfull')
                T = S.Tfull;
                source = 'MAT:Tfull';

            elseif isfield(S, 'Treview')
                T = S.Treview;
                source = 'MAT:Treview';

            else
                f = fieldnames(S);
                T = [];
                source = '';

                for i = 1:numel(f)
                    if istable(S.(f{i}))
                        T = S.(f{i});
                        source = ['MAT:' f{i}];
                        break;
                    end
                end

                if isempty(T)
                    error('MAT review file does not contain a table.');
                end
            end

        otherwise

            T = readtable(reviewFile);
            source = ['TABLE:' ext];
    end

    if ~istable(T)
        error('Loaded review object is not a table.');
    end
end

function name = find_var_name(T, aliases)

    vars = T.Properties.VariableNames;
    normVars = normalize_names(vars);
    normAliases = normalize_names(aliases);

    name = '';

    for i = 1:numel(normAliases)
        idx = find(normVars == normAliases(i), 1, 'first');
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

    % Match Step 6 / model preprocessing.
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

function [centerSampleLocal, sourceName] = get_candidate_local_samples(T, Fs, fs_hdr, first_samp)

    % Priority 1: local sample columns saved by Step 6 full MAT.
    name = find_var_name(T, { ...
        'refinedSampleLocal', ...
        'centerSampleLocal', ...
        'modelCenterSampleLocal'});

    if ~isempty(name)
        centerSampleLocal = round(col_to_double(T.(name)));
        sourceName = name;
        return;
    end

    % Priority 2: absolute FIFF/mbrowse sample columns.
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

    % Priority 3: local time columns in seconds.
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

    % Priority 4: absolute FIFF/mbrowse time columns.
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