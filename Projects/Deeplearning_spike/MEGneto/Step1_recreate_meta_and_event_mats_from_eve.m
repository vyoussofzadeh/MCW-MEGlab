%% ------------------------------------------------------------------------
% Recreate reviewer meta_files + Spike/NoSpike .mat files from .eve files
%
% Input:
%   review_assignments.csv
%   .eve files with spike times
%   matching MEG FIF files
%
% Output:
%   /Processed/Spike-NoSpike_Mat_files/<Rater>/meta_files/*_meta_data.txt
%   /Processed/Spike-NoSpike_Mat_files/<Rater>/*_spike_events.mat
%   /Processed/Spike-NoSpike_Mat_files/<Rater>/*_nospike_events.mat
%
% Events are saved as: trials x time x channels
% -------------------------------------------------------------------------

clc; clear;

%% ======================== USER SETTINGS =================================

% opts.ftRoot = '/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';

opts.ftRoot = '/MEG_data/MEG_Tools/fieldtrip/fieldtrip_2022';

opts.assignCsv = '/MEG_data/AHW_SpikeAnalysis/Review_sheet/review_assignments - Adi.csv';

opts.eveDir = '/MEG_data/AHW_SpikeAnalysis/MEG_data/sss/Spike_times';

opts.outRoot = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files';

opts.megSearchDirs = {
    '/MEG_data/AHW_SpikeAnalysis/MEG_data/sss'
    '/MEG_data/AHW_SpikeAnalysis/MEG_data'
    };

opts.raters = {'Adi','Josh','Manoj','Pradeep'};

% Your later script crops to 306 samples. This window gives 306 samples at 200 Hz:
% -0.500 to +1.025 sec = 100 pre + 205 post + center = 306
opts.targetFs = 200;
opts.epochWinSec = [-0.500 1.025];

opts.keepFirstNChannels = 306;
opts.useSingle = true;

% No-spike generation
opts.makeNoSpike = true;
opts.nNoSpikePerSpike = 1;
opts.excludeAroundSpikeSec = 1.0;
opts.minNoSpikeSeparationSec = 1.0;

opts.overwrite = true;   % set false after first successful recreation
opts.dryRun = false;

rng(20260119);

mneRoot = '/MEG_data/MEG_Tools/mne_matlab_20226/mne-matlab-master/mne-matlab-master/matlab';
addpath(genpath(mneRoot));
% rehash toolboxcache;

fprintf('\nUsing MNE/FIFF reader from:\n%s\n', mneRoot);
which fiff_find_evoked -all
which fiff_setup_read_raw -all

%% ======================== INITIALIZE =====================================

if exist(opts.ftRoot, 'dir')
    addpath(opts.ftRoot);
    ft_defaults;
else
    error('FieldTrip folder not found: %s', opts.ftRoot);
end

if ~isfile(opts.assignCsv)
    error('Missing assignment CSV: %s', opts.assignCsv);
end

if ~exist(opts.outRoot, 'dir')
    mkdir(opts.outRoot);
end

% Create rater folders and meta_files folders
for r = 1:numel(opts.raters)
    raterDir = fullfile(opts.outRoot, opts.raters{r});
    metaDir  = fullfile(raterDir, 'meta_files');

    if ~exist(raterDir, 'dir')
        mkdir(raterDir);
    end

    if ~exist(metaDir, 'dir')
        mkdir(metaDir);
    end
end

Tassign = readtable(opts.assignCsv);

requiredCols = {'FileName','FullPath','Physician1','Physician2'};
for c = 1:numel(requiredCols)
    if ~ismember(requiredCols{c}, Tassign.Properties.VariableNames)
        error('Assignment CSV is missing column: %s', requiredCols{c});
    end
end

fprintf('Loaded assignment CSV: %d rows\n', height(Tassign));

Log = {};

%% ======================== MAIN LOOP ======================================

for i = 1:height(Tassign)

    eveFile = get_table_str(Tassign, i, 'FullPath');

    if isempty(eveFile) || ~isfile(eveFile)
        fname = get_table_str(Tassign, i, 'FileName');
        eveFile = fullfile(opts.eveDir, fname);
    end

    if ~isfile(eveFile)
        warning('Missing .eve file for row %d: %s', i, eveFile);
        continue;
    end

    [~, base, ~] = fileparts(eveFile);

    reviewers = unique({
        get_table_str(Tassign, i, 'Physician1')
        get_table_str(Tassign, i, 'Physician2')
        }, 'stable');

    reviewers = reviewers(~cellfun(@isempty, reviewers));
    reviewers = intersect(reviewers, opts.raters, 'stable');

    if isempty(reviewers)
        warning('No valid reviewers for: %s', base);
        continue;
    end

    fprintf('\n====================================================\n');
    fprintf('Row %d/%d: %s\n', i, height(Tassign), base);
    fprintf('Reviewers: %s\n', strjoin(reviewers, ', '));

    %% -------- Read .eve spike events -------------------------------------

    try
        [sampleRaw, code] = read_eve_file_numeric(eveFile);
    catch ME
        warning('Could not read .eve file:\n%s\n%s', eveFile, ME.message);
        continue;
    end

    if isempty(sampleRaw)
        warning('No events found in: %s', eveFile);
        continue;
    end

    fprintf('Events in .eve: %d\n', numel(sampleRaw));

    %% -------- Find and load MEG file -------------------------------------

    rawFile = find_matching_meg_file(base, opts.megSearchDirs);

    if isempty(rawFile)
        warning('Could not find matching FIF file for base: %s', base);
        continue;
    end

    fprintf('MEG file: %s\n', rawFile);

    try
        [x, chanLabels, Fs, hdrInfo] = load_meg_matrix(rawFile, opts);
    catch ME
        warning('Could not load MEG file:\n%s\n%s', rawFile, ME.message);
        continue;
    end

    nData = size(x, 2);
    fprintf('Loaded MEG: %d channels x %d samples, Fs = %.3f Hz\n', ...
        size(x,1), size(x,2), Fs);

    %% -------- Map .eve samples to loaded/resampled data samples -----------

    sample_ds = map_event_samples_auto(sampleRaw, hdrInfo, Fs, nData);

    nPre  = round(abs(opts.epochWinSec(1)) * Fs);
    nPost = round(abs(opts.epochWinSec(2)) * Fs);
    tsec  = (-nPre:nPost) ./ Fs;

    validEpoch = sample_ds > nPre & sample_ds <= (nData - nPost);

    if ~any(validEpoch)
        warning('No valid epochs after boundary check for %s', base);
        continue;
    end

    spikeCenters = sample_ds(validEpoch);
    spikeCode    = code(validEpoch);
    spikeRawSamp = sampleRaw(validEpoch);

    fprintf('Valid spike epochs: %d / %d\n', numel(spikeCenters), numel(sampleRaw));

    %% -------- Extract spike epochs ---------------------------------------

    [spikeEvents, keptSpike] = extract_epochs(x, spikeCenters, nPre, nPost, opts.useSingle);

    spikeCenters = spikeCenters(keptSpike);
    spikeCode    = spikeCode(keptSpike);
    spikeRawSamp = spikeRawSamp(keptSpike);

    %% -------- Extract random no-spike epochs -----------------------------

    noEvents = [];
    noCenters = [];

    if opts.makeNoSpike
        nNeed = max(1, round(numel(spikeCenters) * opts.nNoSpikePerSpike));

        noCenters = pick_nospike_centers( ...
            nData, nNeed, sample_ds(validEpoch), nPre, nPost, Fs, opts);

        if ~isempty(noCenters)
            [noEvents, keptNo] = extract_epochs(x, noCenters, nPre, nPost, opts.useSingle);
            noCenters = noCenters(keptNo);
        end

        fprintf('No-spike epochs: %d\n', size(noEvents,1));
    end

    %% -------- Save files for each assigned reviewer ----------------------

    for rr = 1:numel(reviewers)

        rater = reviewers{rr};

        raterDir = fullfile(opts.outRoot, rater);
        metaDir  = fullfile(raterDir, 'meta_files');

        outSpike   = fullfile(raterDir, [base '_spike_events.mat']);
        outNoSpike = fullfile(raterDir, [base '_nospike_events.mat']);
        outMeta    = fullfile(metaDir, sprintf('%s_%s_meta_data.txt', base, rater));

        if ~opts.overwrite && isfile(outSpike) && isfile(outNoSpike) && isfile(outMeta)
            fprintf('Already exists, skipping rater %s: %s\n', rater, base);
            continue;
        end

        % Meta file: one row per .eve event.
        % agree=1 only if epoch is valid/extractable.
        agree = double(validEpoch(:));

        Tmeta = table();
        Tmeta.sample      = sampleRaw(:);
        Tmeta.agree       = agree(:);
        Tmeta.code        = code(:);
        Tmeta.time_sec    = double(sample_ds(:) - 1) ./ Fs;
        Tmeta.sample_ds   = sample_ds(:);
        Tmeta.valid_epoch = validEpoch(:);

        if ~opts.dryRun
            writetable(Tmeta, outMeta, 'FileType', 'text', 'Delimiter', '\t');
        end

        % Save spike mat
        events = spikeEvents;
        metas = struct();
        metas.Fs = Fs;
        metas.tsec = tsec;
        metas.winSec = opts.epochWinSec;
        metas.chanLabels = chanLabels;
        metas.base = base;
        metas.rater = rater;
        metas.rawFile = rawFile;
        metas.eveFile = eveFile;
        metas.metaFile = outMeta;
        metas.sample = spikeRawSamp(:);
        metas.sample_ds = spikeCenters(:);
        metas.code = spikeCode(:);
        metas.label = ones(numel(spikeCenters), 1);
        metas.kind = 'spike_from_eve';

        if ~opts.dryRun
            save(outSpike, 'events', 'metas', '-v7.3');
        end

        % Save no-spike mat
        if opts.makeNoSpike && ~isempty(noEvents)
            events = noEvents;
            metas = struct();
            metas.Fs = Fs;
            metas.tsec = tsec;
            metas.winSec = opts.epochWinSec;
            metas.chanLabels = chanLabels;
            metas.base = base;
            metas.rater = rater;
            metas.rawFile = rawFile;
            metas.eveFile = eveFile;
            metas.metaFile = outMeta;
            metas.sample = nan(numel(noCenters), 1);
            metas.sample_ds = noCenters(:);
            metas.code = zeros(numel(noCenters), 1);
            metas.label = zeros(numel(noCenters), 1);
            metas.kind = 'nospike_random';

            if ~opts.dryRun
                save(outNoSpike, 'events', 'metas', '-v7.3');
            end
        end

        fprintf('Saved for %s: %s\n', rater, base);

        Log(end+1,:) = {rater, base, eveFile, rawFile, outMeta, outSpike, outNoSpike, ...
            numel(sampleRaw), sum(validEpoch), size(spikeEvents,1), size(noEvents,1)}; %#ok<SAGROW>
    end
end

%% ======================== SAVE LOG =======================================

if ~isempty(Log)
    Tlog = cell2table(Log, 'VariableNames', ...
        {'Rater','Base','EveFile','RawFile','MetaFile','SpikeMat','NoSpikeMat', ...
         'NEveEvents','NValidMeta','NSpikeSaved','NNoSpikeSaved'});

    disp(Tlog);

    logFile = fullfile(opts.outRoot, ...
        ['recreate_meta_and_mats_log_' datestr(now,'yyyymmdd_HHMMSS') '.csv']);

    if ~opts.dryRun
        writetable(Tlog, logFile);
        fprintf('\nSaved log: %s\n', logFile);
    end
else
    fprintf('\nNo files were generated.\n');
end

%% ======================== LOCAL FUNCTIONS ================================

function s = get_table_str(T, row, colName)
    v = T.(colName);

    if iscell(v)
        s = v{row};
    elseif isstring(v)
        s = v(row);
    elseif iscategorical(v)
        s = char(v(row));
    else
        s = v(row);
    end

    s = char(string(s));
    s = strtrim(s);

    if strcmpi(s, '<missing>') || strcmpi(s, 'missing') || strcmpi(s, 'nan')
        s = '';
    end
end

function [sample, code] = read_eve_file_numeric(eveFile)

    fid = fopen(eveFile, 'r');
    if fid < 0
        error('Cannot open file: %s', eveFile);
    end

    C = {};
    while true
        line = fgetl(fid);
        if ~ischar(line), break; end

        line = strtrim(line);
        if isempty(line), continue; end
        if startsWith(line, '#'), continue; end

        nums = sscanf(line, '%f');
        if isempty(nums), continue; end

        C{end+1,1} = nums(:)'; %#ok<AGROW>
    end

    fclose(fid);

    if isempty(C)
        sample = [];
        code = [];
        return;
    end

    maxLen = max(cellfun(@numel, C));
    M = nan(numel(C), maxLen);

    for i = 1:numel(C)
        M(i,1:numel(C{i})) = C{i};
    end

    sample = M(:,1);

    % MNE-style .eve usually has: sample, previous code, event code
    if size(M,2) >= 3
        code = M(:,3);
    elseif size(M,2) >= 2
        code = M(:,2);
    else
        code = ones(size(sample));
    end

    good = ~isnan(sample);
    sample = sample(good);
    code = code(good);

    code(isnan(code)) = 1;
end

function rawFile = find_matching_meg_file(base, searchDirs)

    rawFile = '';

    exactPatterns = {
        [base '.fif']
        [base '.fif.gz']
        [base '_raw.fif']
        [base '_raw.fif.gz']
        [base '*fif*']
        };

    for d = 1:numel(searchDirs)
        root = searchDirs{d};
        if ~exist(root, 'dir'), continue; end

        for p = 1:numel(exactPatterns)
            L = dir(fullfile(root, '**', exactPatterns{p}));
            if ~isempty(L)
                rawFile = fullfile(L(1).folder, L(1).name);
                return;
            end
        end
    end

    % Fallback: match subject + run number
    m = regexpi(base, '^(?<sub>.+?)_Run0?(?<run>\d+)', 'names', 'once');

    if isempty(m)
        m = regexpi(base, '^(?<sub>.+?)_run0?(?<run>\d+)', 'names', 'once');
    end

    if ~isempty(m)
        rnum = str2double(m.run);
        patterns = {
            sprintf('*%s*Run%02d*fif*', m.sub, rnum)
            sprintf('*%s*run%02d*fif*', m.sub, rnum)
            sprintf('*%s*Run%d*fif*',   m.sub, rnum)
            sprintf('*%s*run%d*fif*',   m.sub, rnum)
            };

        for d = 1:numel(searchDirs)
            root = searchDirs{d};
            if ~exist(root, 'dir'), continue; end

            for p = 1:numel(patterns)
                L = dir(fullfile(root, '**', patterns{p}));
                if ~isempty(L)
                    rawFile = fullfile(L(1).folder, L(1).name);
                    return;
                end
            end
        end
    end
end

function [x, chanLabels, Fs, hdrInfo] = load_meg_matrix(rawFile, opts)

    hdr = ft_read_header(rawFile);

    cfg = [];
    cfg.dataset = rawFile;
    cfg.continuous = 'yes';
    cfg.channel = 'MEG';

    data = ft_preprocessing(cfg);

    Fs0 = data.fsample;

    if ~isempty(opts.targetFs) && abs(Fs0 - opts.targetFs) > 1e-6
        cfg = [];
        cfg.resamplefs = opts.targetFs;
        cfg.detrend = 'no';
        cfg.demean = 'no';
        data = ft_resampledata(cfg, data);
    end

    x = data.trial{1};
    chanLabels = data.label(:);
    Fs = data.fsample;

    if ~isempty(opts.keepFirstNChannels) && size(x,1) > opts.keepFirstNChannels
        x = x(1:opts.keepFirstNChannels, :);
        chanLabels = chanLabels(1:opts.keepFirstNChannels);
    end

    hdrInfo = struct();
    hdrInfo.FsRaw = hdr.Fs;
    hdrInfo.nSamplesRaw = hdr.nSamples;
    hdrInfo.firstSample = get_first_sample(hdr);
end

function firstSample = get_first_sample(hdr)

    firstSample = 0;

    if isfield(hdr, 'FirstSample')
        firstSample = double(hdr.FirstSample);
        return;
    end

    if isfield(hdr, 'firstSample')
        firstSample = double(hdr.firstSample);
        return;
    end

    try
        firstSample = double(hdr.orig.raw.first_samp);
        return;
    catch
    end

    try
        firstSample = double(hdr.orig.first_samp);
        return;
    catch
    end
end

function sample_ds = map_event_samples_auto(sampleRaw, hdrInfo, Fs, nData)

    sampleRaw = double(sampleRaw(:));

    % Option 1: .eve samples are already in loaded-data sample units
    sampleAsData = round(sampleRaw);

    % Option 2: .eve samples are raw FIF sample indices
    eventSec = (sampleRaw - double(hdrInfo.firstSample)) ./ double(hdrInfo.FsRaw);
    sampleAsRaw = round(eventSec .* Fs) + 1;

    validData = sum(sampleAsData >= 1 & sampleAsData <= nData);
    validRaw  = sum(sampleAsRaw  >= 1 & sampleAsRaw  <= nData);

    % Heuristic:
    % If original sample numbers fit directly in the resampled data length,
    % use them directly. Otherwise use raw FIF conversion.
    if max(sampleRaw) <= nData && validData >= validRaw
        sample_ds = sampleAsData;
        fprintf('Event mapping: using .eve samples directly as data samples\n');
    else
        sample_ds = sampleAsRaw;
        fprintf('Event mapping: converting .eve raw samples using hdr Fs/firstSample\n');
    end

    fprintf('Valid by direct mapping: %d | valid by raw conversion: %d\n', ...
        validData, validRaw);
end

function [E, kept] = extract_epochs(x, centers, nPre, nPost, useSingle)

    centers = round(centers(:));
    nChan = size(x, 1);
    nData = size(x, 2);
    nSamp = nPre + nPost + 1;

    valid = centers > nPre & centers <= (nData - nPost);
    centers = centers(valid);

    if useSingle
        E = zeros(numel(centers), nSamp, nChan, 'single');
    else
        E = zeros(numel(centers), nSamp, nChan);
    end

    for i = 1:numel(centers)
        idx = (centers(i)-nPre):(centers(i)+nPost);
        seg = x(:, idx)';  % time x channel

        if useSingle
            E(i,:,:) = single(seg);
        else
            E(i,:,:) = seg;
        end
    end

    kept = valid;
end

function noCenters = pick_nospike_centers(nData, nNeed, spikeCenters, nPre, nPost, Fs, opts)

    valid = true(nData, 1);

    valid(1:nPre) = false;
    valid((nData-nPost+1):end) = false;

    excludeSamp = round(opts.excludeAroundSpikeSec * Fs);
    spikeCenters = round(spikeCenters(:));
    spikeCenters = spikeCenters(spikeCenters >= 1 & spikeCenters <= nData);

    for i = 1:numel(spikeCenters)
        lo = max(1, spikeCenters(i) - excludeSamp);
        hi = min(nData, spikeCenters(i) + excludeSamp);
        valid(lo:hi) = false;
    end

    candidates = find(valid);

    if isempty(candidates)
        noCenters = [];
        return;
    end

    candidates = candidates(randperm(numel(candidates)));

    minSep = round(opts.minNoSpikeSeparationSec * Fs);

    noCenters = [];

    for i = 1:numel(candidates)
        c = candidates(i);

        if isempty(noCenters) || all(abs(c - noCenters) > minSep)
            noCenters(end+1,1) = c; %#ok<AGROW>
        end

        if numel(noCenters) >= nNeed
            break;
        end
    end

    if numel(noCenters) < nNeed
        warning('Requested %d no-spike epochs, found %d.', nNeed, numel(noCenters));
    end
end