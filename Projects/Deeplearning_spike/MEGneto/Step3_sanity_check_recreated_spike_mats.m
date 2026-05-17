%% ------------------------------------------------------------------------
% Sanity check recreated spike/no-spike .mat files and meta files
% -------------------------------------------------------------------------

clc; clear;

root = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files';
raters = {'Adi','Josh','Manoj','Pradeep'};

expectedFs = 200;
expectedNTime = 306;
expectedNChan = 306;

Rows = {};
Errors = {};

for r = 1:numel(raters)

    rater = raters{r};
    raterDir = fullfile(root, rater);
    metaDir  = fullfile(raterDir, 'meta_files');

    fprintf('\n==============================\n');
    fprintf('Rater: %s\n', rater);
    fprintf('==============================\n');

    spkFiles = dir(fullfile(raterDir, '*_spike_events.mat'));
    nosFiles = dir(fullfile(raterDir, '*_nospike_events.mat'));
    metaFiles = dir(fullfile(metaDir, sprintf('*_%s_meta_data.txt', rater)));

    fprintf('Spike mats:    %d\n', numel(spkFiles));
    fprintf('No-spike mats: %d\n', numel(nosFiles));
    fprintf('Meta files:    %d\n', numel(metaFiles));

    if isempty(spkFiles)
        Errors(end+1,:) = {rater, '', 'No spike mat files found'}; %#ok<SAGROW>
        continue;
    end

    for k = 1:numel(spkFiles)

        spkFn = fullfile(spkFiles(k).folder, spkFiles(k).name);
        base = erase(spkFiles(k).name, '_spike_events.mat');

        nosFn  = fullfile(raterDir, [base '_nospike_events.mat']);
        metaFn = fullfile(metaDir, sprintf('%s_%s_meta_data.txt', base, rater));

        hasNoSpike = isfile(nosFn);
        hasMeta    = isfile(metaFn);

        nSpike = NaN;
        nNoSpike = NaN;
        nAgree = NaN;
        nMeta = NaN;
        nTime = NaN;
        nChan = NaN;
        Fs = NaN;
        hasNaN = false;
        hasInf = false;
        medAbs = NaN;

        status = "OK";

        %% ---- Check required companion files ----
        if ~hasMeta
            Errors(end+1,:) = {rater, base, 'Missing meta file'}; %#ok<SAGROW>
            status = "ERROR";
        end

        if ~hasNoSpike
            Errors(end+1,:) = {rater, base, 'Missing no-spike mat file'}; %#ok<SAGROW>
            status = "ERROR";
        end

        %% ---- Load spike mat ----
        try
            [E, M] = load_event_mat(spkFn);

            sz = size(E);
            if numel(sz) < 3
                Errors(end+1,:) = {rater, base, 'Spike events is not 3D'}; %#ok<SAGROW>
                status = "ERROR";
            else
                nSpike = sz(1);
                nTime  = sz(2);
                nChan  = sz(3);
            end

            hasNaN = any(isnan(E(:)));
            hasInf = any(isinf(E(:)));
            medAbs = median(abs(double(E(:))), 'omitnan');

            if isstruct(M) && isfield(M, 'Fs')
                Fs = double(M.Fs);
                if numel(Fs) > 1
                    Fs = Fs(1);
                end
            end

            if hasNaN
                Errors(end+1,:) = {rater, base, 'Spike events contain NaN'}; %#ok<SAGROW>
                status = "ERROR";
            end

            if hasInf
                Errors(end+1,:) = {rater, base, 'Spike events contain Inf'}; %#ok<SAGROW>
                status = "ERROR";
            end

            if nTime ~= expectedNTime
                Errors(end+1,:) = {rater, base, sprintf('Unexpected nTime: %d', nTime)}; %#ok<SAGROW>
                status = "ERROR";
            end

            if nChan ~= expectedNChan
                Errors(end+1,:) = {rater, base, sprintf('Unexpected nChan: %d', nChan)}; %#ok<SAGROW>
                status = "ERROR";
            end

            if ~isnan(Fs) && abs(Fs - expectedFs) > 1e-6
                Errors(end+1,:) = {rater, base, sprintf('Unexpected Fs: %.3f', Fs)}; %#ok<SAGROW>
                status = "ERROR";
            end

            if medAbs == 0 || isnan(medAbs)
                Errors(end+1,:) = {rater, base, 'Median absolute signal is zero or NaN'}; %#ok<SAGROW>
                status = "WARNING";
            end

        catch ME
            Errors(end+1,:) = {rater, base, ['Could not load spike mat: ' ME.message]}; %#ok<SAGROW>
            status = "ERROR";
        end

        %% ---- Load no-spike mat ----
        if hasNoSpike
            try
                [En, Mn] = load_event_mat(nosFn); %#ok<ASGLU>
                nNoSpike = size(En,1);

                if size(En,2) ~= expectedNTime || size(En,3) ~= expectedNChan
                    Errors(end+1,:) = {rater, base, ...
                        sprintf('No-spike size mismatch: [%d %d %d]', size(En,1), size(En,2), size(En,3))}; %#ok<SAGROW>
                    status = "ERROR";
                end

                if any(isnan(En(:))) || any(isinf(En(:)))
                    Errors(end+1,:) = {rater, base, 'No-spike events contain NaN or Inf'}; %#ok<SAGROW>
                    status = "ERROR";
                end

            catch ME
                Errors(end+1,:) = {rater, base, ['Could not load no-spike mat: ' ME.message]}; %#ok<SAGROW>
                status = "ERROR";
            end
        end

        %% ---- Read meta file ----
        if hasMeta
            try
                Tm = read_meta_flexible(metaFn);
                nMeta = height(Tm);

                required = {'sample','agree','code'};
                if ~all(ismember(required, Tm.Properties.VariableNames))
                    Errors(end+1,:) = {rater, base, 'Meta missing sample/agree/code'}; %#ok<SAGROW>
                    status = "ERROR";
                else
                    nAgree = sum(round(Tm.agree) == 1);

                    if nSpike ~= nAgree
                        Errors(end+1,:) = {rater, base, ...
                            sprintf('Spike/meta mismatch: nSpike=%d, meta agree==1=%d', nSpike, nAgree)}; %#ok<SAGROW>
                        status = "ERROR";
                    end
                end

            catch ME
                Errors(end+1,:) = {rater, base, ['Could not read meta: ' ME.message]}; %#ok<SAGROW>
                status = "ERROR";
            end
        end

        Rows(end+1,:) = {rater, base, status, nSpike, nNoSpike, nMeta, nAgree, ...
            nTime, nChan, Fs, hasMeta, hasNoSpike, hasNaN, hasInf, medAbs, ...
            spkFn, metaFn, nosFn}; %#ok<SAGROW>

        clear E En M Mn Tm
    end
end

%% ---- Summary table ----

Tcheck = cell2table(Rows, 'VariableNames', ...
    {'Rater','Base','Status','NSpike','NNoSpike','NMeta','NAgree', ...
     'NTime','NChan','Fs','HasMeta','HasNoSpike','HasNaN','HasInf','MedianAbs', ...
     'SpikeFile','MetaFile','NoSpikeFile'});

disp(Tcheck(:, {'Rater','Base','Status','NSpike','NNoSpike','NMeta','NAgree','NTime','NChan','Fs'}));

fprintf('\nTotal files checked: %d\n', height(Tcheck));
fprintf('OK:      %d\n', sum(Tcheck.Status == "OK"));
fprintf('ERROR:   %d\n', sum(Tcheck.Status == "ERROR"));
fprintf('WARNING: %d\n', sum(Tcheck.Status == "WARNING"));

%% ---- Error table ----

if ~isempty(Errors)
    Terr = cell2table(Errors, 'VariableNames', {'Rater','Base','Issue'});
    fprintf('\nProblems found:\n');
    disp(Terr);

    outErr = fullfile(root, ['sanity_check_errors_' datestr(now,'yyyymmdd_HHMMSS') '.csv']);
    writetable(Terr, outErr);
    fprintf('Saved error report: %s\n', outErr);
else
    fprintf('\nNo errors found.\n');
end

outCheck = fullfile(root, ['sanity_check_summary_' datestr(now,'yyyymmdd_HHMMSS') '.csv']);
writetable(Tcheck, outCheck);
fprintf('Saved summary report: %s\n', outCheck);

%% ---- Plot one example for visual check ----

idxOK = find(Tcheck.Status == "OK", 1, 'first');

if ~isempty(idxOK)
    [E, M] = load_event_mat(Tcheck.SpikeFile{idxOK});

    trial = 1;

    if isstruct(M) && isfield(M, 'tsec')
        t = M.tsec;
    else
        t = (0:size(E,2)-1) ./ expectedFs;
    end

    figure;
    imagesc(t, 1:size(E,3), squeeze(E(trial,:,:))');
    axis xy tight;
    xlabel('Time (s)');
    ylabel('Channel');
    title(sprintf('Example spike epoch: %s | %s | trial %d', ...
        Tcheck.Rater{idxOK}, Tcheck.Base{idxOK}, trial), 'Interpreter', 'none');
    colorbar;

    figure;
    plot(t, squeeze(E(trial,:,1)));
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Example spike epoch, channel 1');
    grid on;
end

%% ======================= helper functions ===============================

function [E, M] = load_event_mat(fn)

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
        opts = detectImportOptions(fn, 'FileType', 'text');
        T = readtable(fn, opts);
    catch
        T = readtable(fn, 'FileType', 'text', 'Delimiter', '\t');
    end

    % Normalize variable names, just in case
    oldNames = T.Properties.VariableNames;
    newNames = lower(regexprep(oldNames, '[^a-zA-Z0-9]', ''));

    for i = 1:numel(newNames)
        switch newNames{i}
            case {'sample','sampleidx','sampleindex','eventsample'}
                T.Properties.VariableNames{i} = 'sample';
            case {'agree','accepted','accept','include','included','isspike','spike','label'}
                T.Properties.VariableNames{i} = 'agree';
            case {'code','eventcode','eventid','trigger','value'}
                T.Properties.VariableNames{i} = 'code';
        end
    end

    if ismember('agree', T.Properties.VariableNames)
        T.agree = double(T.agree);
    end

    if ismember('sample', T.Properties.VariableNames)
        T.sample = double(T.sample);
    end

    if ismember('code', T.Properties.VariableNames)
        T.code = double(T.code);
    end
end