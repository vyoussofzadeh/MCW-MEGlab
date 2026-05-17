%% ------------------------------------------------------------------------
% Recreate per-reviewer Spike/NoSpike event .mat files from reviewer meta txt
%
% Expected output:
%   /MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files/Adi/
%       <base>_spike_events.mat
%       <base>_nospike_events.mat
%       meta_files/<base>_Adi_meta_data.txt
%
% The later consensus script expects spike epochs as trials x time x channels.
% This script saves:
%       events : N x time x channels
%       metas  : struct with Fs, tsec, chanLabels, samples, etc.
% -------------------------------------------------------------------------

clc; clear;

%% ======================== USER SETTINGS =================================

opts.ftRoot = '/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';

opts.outRoot = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files';

opts.megSearchDirs = {
    '/MEG_data/AHW_SpikeAnalysis/MEG_data/sss'
    '/MEG_data/AHW_SpikeAnalysis/MEG_data'
    };

opts.raters = {'Adi','Josh','Manoj','Pradeep'};

% Output epoch sampling rate
opts.targetFs = 200;

% Epoch window around spike center.
% At 200 Hz, [-0.500 1.025] gives 306 samples:
% 100 pre + 205 post + center = 306
opts.epochWinSec = [-0.500 1.025];

% If meta sample values came from .eve/raw FIF sample numbers, use 'raw'.
% If they are already in the 200 Hz exported data sample index, use 'target'.
% If your meta file has time_sec/time column, 'auto' will prefer time.
opts.metaSampleRate = 'auto';   % 'auto', 'raw', 'target', or 'time'

% No-spike generation
opts.makeNoSpike = true;
opts.nNoSpikePerSpike = 1;          % usually 1 negative per positive
opts.excludeAroundEventsSec = 1.0;  % avoid all candidate events by +/- this
opts.minNoSpikeSeparationSec = 1.0;

% File behavior
opts.overwriteMat = false;
opts.rewriteMetaAsStandardTxt = true;   % keeps sample/agree/code readable
opts.backupOldMeta = true;
opts.dryRun = false;

% Data handling
opts.useSingle = true;
opts.keepFirstNChannels = 306;

rng(20260119);

%% ======================== INITIALIZE =====================================

if ~isempty(opts.ftRoot) && exist(opts.ftRoot, 'dir')
    addpath(opts.ftRoot);
    ft_defaults;
end

assert(exist('ft_preprocessing','file') == 2, ...
    'FieldTrip not found. Check opts.ftRoot.');

Log = {};

%% ======================== MAIN LOOP ======================================

for rr = 1:numel(opts.raters)

    rater = opts.raters{rr};
    raterDir = fullfile(opts.outRoot, rater);
    metaDir  = fullfile(raterDir, 'meta_files');

    if ~exist(metaDir, 'dir')
        warning('Missing meta folder: %s', metaDir);
        continue;
    end

    metaFiles = dir(fullfile(metaDir, sprintf('*_%s_meta_data.txt', rater)));

    fprintf('\n==============================\n');
    fprintf('Rater: %s | meta files: %d\n', rater, numel(metaFiles));
    fprintf('==============================\n');

    for mm = 1:numel(metaFiles)

        metaFile = fullfile(metaFiles(mm).folder, metaFiles(mm).name);

        suffix = sprintf('_%s_meta_data.txt', rater);
        base = erase(metaFiles(mm).name, suffix);

        outSpike   = fullfile(raterDir, [base '_spike_events.mat']);
        outNoSpike = fullfile(raterDir, [base '_nospike_events.mat']);

        needSpike   = opts.overwriteMat || ~isfile(outSpike);
        needNoSpike = opts.makeNoSpike && (opts.overwriteMat || ~isfile(outNoSpike));

        if ~needSpike && ~needNoSpike
            fprintf('Already exists, skipping: %s\n', base);
            continue;
        end

        fprintf('\n--- %s | %s ---\n', rater, base);

        try
            Tmeta = read_meta_table_flexible(metaFile);
        catch ME
            warning('Could not read meta file:\n%s\n%s', metaFile, ME.message);
            continue;
        end

        if isempty(Tmeta) || ~all(ismember({'sample','agree','code'}, Tmeta.Properties.VariableNames))
            warning('Meta file lacks required sample/agree/code columns: %s', metaFile);
            continue;
        end

        posRows = find(round(Tmeta.agree) == 1);
        if isempty(posRows)
            warning('No agree==1 rows in meta: %s', metaFile);
            continue;
        end

        rawFile = find_matching_meg_file(base, opts.megSearchDirs);
        if isempty(rawFile)
            warning('Could not find matching MEG FIF for base: %s', base);
            continue;
        end

        fprintf('MEG file: %s\n', rawFile);

        try
            [x, chanLabels, Fs, hdrInfo] = load_meg_as_matrix(rawFile, opts);
        catch ME
            warning('Could not load MEG file:\n%s\n%s', rawFile, ME.message);
            continue;
        end

        nTimeData = size(x, 2);

        % Map meta sample or time values to sample indices in the loaded/resampled data
        sample_ds_all = map_meta_samples_to_data_indices(Tmeta, hdrInfo, Fs, nTimeData, opts);

        nPre  = round(abs(opts.epochWinSec(1)) * Fs);
        nPost = round(abs(opts.epochWinSec(2)) * Fs);
        tsec  = (-nPre:nPost) ./ Fs;

        validEpoch = sample_ds_all > nPre & sample_ds_all <= (nTimeData - nPost);

        % Approved spikes only, and only if extractable
        posRowsValid = posRows(validEpoch(posRows));
        spikeCenters = sample_ds_all(posRowsValid);

        if isempty(spikeCenters)
            warning('No valid spike epochs after boundary check: %s', base);
            continue;
        end

        % Rewrite standardized meta so approved rows match saved spike epochs
        if opts.rewriteMetaAsStandardTxt && ~opts.dryRun
            Tstandard = Tmeta;
            Tstandard.sample_ds = sample_ds_all;
            Tstandard.valid_epoch = validEpoch;

            % Approved rows outside epoch boundary are marked as 0 to avoid
            % mismatch between meta agree==1 count and saved event count.
            badPos = round(Tstandard.agree) == 1 & ~Tstandard.valid_epoch;
            Tstandard.agree(badPos) = 0;

            backup_and_write_meta(metaFile, Tstandard, opts);
        end

        %% ---------------- Save spike epochs -------------------------------
        nSpikeSaved = 0;

        if needSpike
            [events, keptMask] = extract_epochs_from_matrix(x, spikeCenters, nPre, nPost, opts.useSingle);

            keptRows = posRowsValid(keptMask);
            keptCenters = spikeCenters(keptMask);

            metas = struct();
            metas.Fs = Fs;
            metas.tsec = tsec;
            metas.winSec = opts.epochWinSec;
            metas.chanLabels = chanLabels;
            metas.base = base;
            metas.rater = rater;
            metas.rawFile = rawFile;
            metas.metaFile = metaFile;
            metas.sample = Tmeta.sample(keptRows);
            metas.sample_ds = keptCenters;
            metas.code = Tmeta.code(keptRows);
            metas.label = ones(numel(keptCenters), 1);
            metas.kind = 'spike';

            nSpikeSaved = size(events, 1);

            fprintf('Spike epochs: %d\n', nSpikeSaved);

            if ~opts.dryRun
                save(outSpike, 'events', 'metas', '-v7.3');
                fprintf('Saved: %s\n', outSpike);
            end
        else
            fprintf('Spike mat already exists: %s\n', outSpike);
        end

        %% ---------------- Save no-spike epochs ----------------------------
        nNoSaved = 0;

        if needNoSpike
            nNeed = max(1, round(numel(spikeCenters) * opts.nNoSpikePerSpike));

            avoidCenters = sample_ds_all(validEpoch);
            noCenters = pick_nospike_centers( ...
                nTimeData, nNeed, avoidCenters, nPre, nPost, Fs, opts);

            if isempty(noCenters)
                warning('Could not find valid no-spike centers for %s', base);
            else
                [events, keptMask] = extract_epochs_from_matrix(x, noCenters, nPre, nPost, opts.useSingle);
                noCenters = noCenters(keptMask);

                metas = struct();
                metas.Fs = Fs;
                metas.tsec = tsec;
                metas.winSec = opts.epochWinSec;
                metas.chanLabels = chanLabels;
                metas.base = base;
                metas.rater = rater;
                metas.rawFile = rawFile;
                metas.metaFile = metaFile;
                metas.sample = nan(numel(noCenters), 1);
                metas.sample_ds = noCenters(:);
                metas.code = zeros(numel(noCenters), 1);
                metas.label = zeros(numel(noCenters), 1);
                metas.kind = 'nospike_random';

                nNoSaved = size(events, 1);

                fprintf('No-spike epochs: %d\n', nNoSaved);

                if ~opts.dryRun
                    save(outNoSpike, 'events', 'metas', '-v7.3');
                    fprintf('Saved: %s\n', outNoSpike);
                end
            end
        elseif opts.makeNoSpike
            fprintf('No-spike mat already exists: %s\n', outNoSpike);
        end

        Log(end+1,:) = {rater, base, rawFile, metaFile, outSpike, outNoSpike, ...
            numel(posRows), numel(posRowsValid), nSpikeSaved, nNoSaved}; %#ok<SAGROW>
    end
end

%% ======================== SAVE LOG =======================================

if ~isempty(Log)
    Tlog = cell2table(Log, 'VariableNames', ...
        {'Rater','Base','RawFile','MetaFile','SpikeMat','NoSpikeMat', ...
         'NAgreeInMeta','NValidAgree','NSpikeSaved','NNoSpikeSaved'});

    disp(Tlog);

    logFile = fullfile(opts.outRoot, ...
        ['recreate_spike_nospike_log_' datestr(now,'yyyymmdd_HHMMSS') '.csv']);

    if ~opts.dryRun
        writetable(Tlog, logFile);
        fprintf('\nSaved log: %s\n', logFile);
    end
else
    fprintf('\nNo files were generated.\n');
end

%% ======================== LOCAL FUNCTIONS ================================

function T = read_meta_table_flexible(fn)

    % Try normal table read first
    try
        optsImport = detectImportOptions(fn, 'FileType', 'text');
        T0 = readtable(fn, optsImport);
    catch
        T0 = readtable(fn, 'FileType', 'text');
    end

    if isempty(T0)
        T = table();
        return;
    end

    rawNames = T0.Properties.VariableNames;
    names = lower(regexprep(rawNames, '[^a-zA-Z0-9]', ''));

    sampleIdx = find_col(names, {'sample','sampleidx','sampleindex','eventsample','latency'});
    agreeIdx  = find_col(names, {'agree','accepted','accept','include','included','isspike','spike','label','decision'});
    codeIdx   = find_col(names, {'code','eventcode','eventid','trigger','value'});
    timeIdx   = find_col(names, {'time','timesec','timeinsec','seconds','sec'});

    % If column names are not useful, assume common numeric format.
    % Minimum expected: sample, agree, code.
    if isempty(sampleIdx) || isempty(agreeIdx)
        if width(T0) >= 3
            sampleIdx = 1;
            agreeIdx  = 2;
            codeIdx   = 3;
        else
            error('Cannot infer sample/agree/code columns from %s', fn);
        end
    end

    sample = col_to_double(T0{:, sampleIdx});
    agree  = col_to_agree(T0{:, agreeIdx});

    if isempty(codeIdx)
        code = ones(height(T0), 1);
    else
        code = col_to_double(T0{:, codeIdx});
    end

    if isempty(timeIdx)
        time_sec = nan(height(T0), 1);
    else
        time_sec = col_to_double(T0{:, timeIdx});
    end

    T = table(sample(:), agree(:), code(:), time_sec(:), ...
        'VariableNames', {'sample','agree','code','time_sec'});

    good = ~isnan(T.sample) & ~isnan(T.agree);
    T = T(good, :);
end

function idx = find_col(names, keys)
    idx = [];
    keys = lower(regexprep(keys, '[^a-zA-Z0-9]', ''));
    for k = 1:numel(keys)
        ii = find(strcmp(names, keys{k}), 1, 'first');
        if ~isempty(ii)
            idx = ii;
            return;
        end
    end
    for k = 1:numel(keys)
        ii = find(contains(names, keys{k}), 1, 'first');
        if ~isempty(ii)
            idx = ii;
            return;
        end
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
end

function x = col_to_agree(v)
    if isnumeric(v)
        x = double(v);
        return;
    end

    s = lower(strtrim(string(v)));
    x = nan(numel(s), 1);

    yesVals = ["1","yes","y","true","t","spike","agree","accepted","accept","include","included"];
    noVals  = ["0","no","n","false","f","nospike","reject","rejected","exclude","excluded"];

    x(ismember(s, yesVals)) = 1;
    x(ismember(s, noVals))  = 0;

    numericVals = str2double(s);
    useNum = isnan(x) & ~isnan(numericVals);
    x(useNum) = numericVals(useNum);
end

function rawFile = find_matching_meg_file(base, searchDirs)

    rawFile = '';

    exactNames = {
        [base '.fif']
        [base '.fif.gz']
        [base '_raw.fif']
        [base '_raw.fif.gz']
        };

    for d = 1:numel(searchDirs)
        root = searchDirs{d};
        if ~exist(root, 'dir'), continue; end

        for e = 1:numel(exactNames)
            L = dir(fullfile(root, '**', exactNames{e}));
            if ~isempty(L)
                rawFile = fullfile(L(1).folder, L(1).name);
                return;
            end
        end
    end

    % Flexible fallback using subject and RunXX
    m = regexpi(base, '^(?<sub>.+?)_Run0?(?<run>\d+)', 'names', 'once');

    if ~isempty(m)
        pat = sprintf('*%s*Run%02d*fif*', m.sub, str2double(m.run));

        for d = 1:numel(searchDirs)
            root = searchDirs{d};
            if ~exist(root, 'dir'), continue; end

            L = dir(fullfile(root, '**', pat));

            if isempty(L)
                pat2 = sprintf('*%s*run%02d*fif*', m.sub, str2double(m.run));
                L = dir(fullfile(root, '**', pat2));
            end

            if ~isempty(L)
                % Prefer files that contain the full base name if possible
                names = {L.name};
                hit = find(contains(names, base), 1, 'first');
                if isempty(hit), hit = 1; end
                rawFile = fullfile(L(hit).folder, L(hit).name);
                return;
            end
        end
    end

    % Last fallback: base as wildcard
    for d = 1:numel(searchDirs)
        root = searchDirs{d};
        if ~exist(root, 'dir'), continue; end

        L = dir(fullfile(root, '**', [base '*fif*']));
        if ~isempty(L)
            rawFile = fullfile(L(1).folder, L(1).name);
            return;
        end
    end
end

function [x, chanLabels, Fs, hdrInfo] = load_meg_as_matrix(rawFile, opts)

    hdr = ft_read_header(rawFile);

    cfg = [];
    cfg.dataset = rawFile;
    cfg.continuous = 'yes';
    cfg.channel = 'MEG';

    data = ft_preprocessing(cfg);

    Fs0 = data.fsample;

    if ~isempty(opts.targetFs) && opts.targetFs > 0 && abs(Fs0 - opts.targetFs) > 1e-6
        cfg = [];
        cfg.resamplefs = opts.targetFs;
        cfg.detrend = 'no';
        cfg.demean = 'no';
        data = ft_resampledata(cfg, data);
    end

    x = data.trial{1};
    chanLabels = data.label(:);
    Fs = data.fsample;

    if ~isempty(opts.keepFirstNChannels) && numel(chanLabels) > opts.keepFirstNChannels
        keep = 1:opts.keepFirstNChannels;
        x = x(keep, :);
        chanLabels = chanLabels(keep);
    end

    hdrInfo = struct();
    hdrInfo.FsRaw = hdr.Fs;
    hdrInfo.nSamplesRaw = hdr.nSamples;
    hdrInfo.firstSample = get_first_sample_from_header(hdr);
end

function firstSample = get_first_sample_from_header(hdr)
    firstSample = 0;

    candidates = {
        {'FirstSample'}
        {'firstSample'}
        {'orig','first_samp'}
        {'orig','raw','first_samp'}
        {'orig','info','sfreq'}
        };

    for c = 1:numel(candidates)
        path = candidates{c};
        try
            val = hdr;
            for p = 1:numel(path)
                val = val.(path{p});
            end

            if isnumeric(val) && isscalar(val)
                % Do not accidentally use sfreq as first sample
                if val < 1e8 && val ~= hdr.Fs
                    firstSample = double(val);
                    return;
                end
            end
        catch
        end
    end
end

function sample_ds = map_meta_samples_to_data_indices(Tmeta, hdrInfo, Fs, nTimeData, opts)

    hasTime = ismember('time_sec', Tmeta.Properties.VariableNames) && ...
              any(~isnan(Tmeta.time_sec));

    mode = lower(opts.metaSampleRate);

    if strcmp(mode, 'auto')
        if hasTime
            mode = 'time';
        else
            mode = 'raw';
        end
    end

    switch mode
        case 'time'
            if ~hasTime
                error('metaSampleRate is time, but no usable time_sec column was found.');
            end
            sample_ds = round(Tmeta.time_sec .* Fs) + 1;

        case 'target'
            sample_ds = round(Tmeta.sample);

        case 'raw'
            eventSec = (double(Tmeta.sample) - double(hdrInfo.firstSample)) ./ double(hdrInfo.FsRaw);
            sample_ds = round(eventSec .* Fs) + 1;

        otherwise
            error('Unknown opts.metaSampleRate: %s', opts.metaSampleRate);
    end

    sample_ds = double(sample_ds(:));

    % Guard against impossible mappings
    if all(sample_ds < 1 | sample_ds > nTimeData)
        warning(['All mapped event samples are outside the data range. ', ...
                 'Try changing opts.metaSampleRate to ''target'' or ''raw''.']);
    end
end

function [E, keptMask] = extract_epochs_from_matrix(x, centers, nPre, nPost, useSingle)

    centers = round(centers(:));
    nChan = size(x, 1);
    nTimeData = size(x, 2);
    nSamp = nPre + nPost + 1;

    valid = centers > nPre & centers <= (nTimeData - nPost);
    centersValid = centers(valid);

    if useSingle
        E = zeros(numel(centersValid), nSamp, nChan, 'single');
    else
        E = zeros(numel(centersValid), nSamp, nChan);
    end

    for i = 1:numel(centersValid)
        idx = (centersValid(i)-nPre):(centersValid(i)+nPost);
        seg = x(:, idx)';  % time x channels
        if useSingle
            E(i,:,:) = single(seg);
        else
            E(i,:,:) = seg;
        end
    end

    keptMask = valid;
end

function noCenters = pick_nospike_centers(nTimeData, nNeed, avoidCenters, nPre, nPost, Fs, opts)

    valid = true(nTimeData, 1);

    valid(1:nPre) = false;
    valid((nTimeData-nPost+1):end) = false;

    avoidRadius = round(opts.excludeAroundEventsSec * Fs);
    avoidCenters = round(avoidCenters(:));
    avoidCenters = avoidCenters(avoidCenters >= 1 & avoidCenters <= nTimeData);

    for i = 1:numel(avoidCenters)
        lo = max(1, avoidCenters(i) - avoidRadius);
        hi = min(nTimeData, avoidCenters(i) + avoidRadius);
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

function backup_and_write_meta(metaFile, Tstandard, opts)

    if opts.backupOldMeta
        bakFile = [metaFile '.bak_' datestr(now,'yyyymmdd_HHMMSS')];
        if ~isfile(bakFile)
            copyfile(metaFile, bakFile);
        end
    end

    % Keep only simple columns that the later read_meta_txt can parse.
    keepNames = {'sample','agree','code','time_sec','sample_ds','valid_epoch'};
    keepNames = keepNames(ismember(keepNames, Tstandard.Properties.VariableNames));

    Tout = Tstandard(:, keepNames);

    writetable(Tout, metaFile, 'FileType', 'text', 'Delimiter', '\t');
end