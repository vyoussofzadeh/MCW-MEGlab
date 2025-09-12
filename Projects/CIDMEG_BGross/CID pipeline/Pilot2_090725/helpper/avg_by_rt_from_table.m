function results = avg_by_rt_from_table(T, subjectDir, varargin)
% Average Brainstorm trials per run using behavioral table T.
% T must have columns: run, trial, rt_category
% Directory layout: subjectDir/*_Run_%02d/data_<run>_trialNNN.mat
%
% Usage:
% results = avg_by_rt_from_table(T, '/path/to/mcwa086_v1', ...
%                 'IncludeTimeout', false, 'DoContrast', true);

% ---- Options ----
p = inputParser;
addParameter(p,'IncludeTimeout',false,@islogical);
addParameter(p,'DoContrast',true,@islogical);
addParameter(p,'Runs',[],@(x)isnumeric(x) || isempty(x));
parse(p,varargin{:});
opt = p.Results;

% rt_category codes
FAST = 1; SLOW = 2; TIMEOUT = 3;

if isempty(opt.Runs)
    runs = unique(T.run(:))';
else
    runs = opt.Runs(:)';
end

runs = runs+1;

results = struct('run',{},'nFast',{},'nSlow',{},'nTimeout',{}, ...
    'avgFastFile',{},'avgSlowFile',{},'contrastFile',{});

for runNo = runs
    % Locate run folder: *_Run_%02d
    pat = sprintf('*_Run_%02d', runNo);
    d = dir(fullfile(subjectDir, pat));
    if isempty(d)
        fprintf('RUN %d: folder not found (%s)\n', runNo, fullfile(subjectDir, pat));
        continue;
    end
    runDir = fullfile(d(1).folder, d(1).name);

    % Subset rows for this run
    R = T(T.run == runNo, :);
    if isempty(R)
        fprintf('RUN %d: no rows in table\n', runNo);
        continue;
    end

    % Helper to build filename
    % mkFile = @(trial) fullfile(runDir, sprintf('data_%d_trial%03d.mat', runNo, trial));
    mkFile = @(trial) fullfile(runDir, sprintf('data_2_trial%03d.mat', trial));


    % Collect trials by category
    trFast   = R.trial(R.rt_category == FAST);
    trSlow   = R.trial(R.rt_category == SLOW);
    trTO     = R.trial(R.rt_category == TIMEOUT);

    % Average helpers
    avgFastFile = '';
    avgSlowFile = '';
    contrastFile = '';

    if ~isempty(trFast)
        avgFastFile = do_average(runDir, runNo, trFast, mkFile, 'RTfast_fromT');
    end
    if ~isempty(trSlow)
        avgSlowFile = do_average(runDir, runNo, trSlow, mkFile, 'RTslow_fromT');
    end
    if opt.IncludeTimeout && ~isempty(trTO)
        do_average(runDir, runNo, trTO, mkFile, 'RTtimeout_fromT'); %#ok<*NASGU>
    end

    % Optional contrast (Fast - Slow) if both exist
    if opt.DoContrast && ~isempty(avgFastFile) && ~isempty(avgSlowFile)
        Sfast = load(avgFastFile);
        Sslow = load(avgSlowFile);
        % Ensure same #channels and samples
        [Sfast, Sslow] = align_FS(Sfast, Sslow);
        Sdiff = Sfast;
        Sdiff.F = Sfast.F - Sslow.F;
        Sdiff.Comment = sprintf('Contrast |Run=%d| |RT=fast-minus-slow|', runNo);
        if isfield(Sdiff,'History') && iscell(Sdiff.History)
            Sdiff.History(end+1,1:3) = {datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
                'avg_by_rt_from_table', 'Fast - Slow'};
        end
        % contrastFile = fullfile(runDir, sprintf('data_%d_contrast_RTfast_minus_RTslow_fromT.mat', runNo));
        contrastFile = fullfile(runDir, sprintf('data_2_contrast_run%02d_RTfast_minus_RTslow_fromT.mat', runNo));
        save(contrastFile, '-struct', 'Sdiff');
        fprintf('RUN %d: wrote %s\n', runNo, contrastFile);
    end

    results(end+1) = struct('run',runNo, ...
        'nFast',numel(trFast), 'nSlow',numel(trSlow), 'nTimeout',numel(trTO), ...
        'avgFastFile',avgFastFile, 'avgSlowFile',avgSlowFile, ...
        'contrastFile',contrastFile);
end
end

% --------- helpers ---------

function outFile = do_average(runDir, runNo, trials, mkFile, tag)
% Load listed trial files, average F, write Brainstorm-compatible MAT
Flist = {};
Slast = [];
missing = 0;

for t = trials(:)'
    fn = mkFile(t);
    if ~isfile(fn)
        fprintf('  MISSING: %s\n', fn);
        missing = missing + 1;
        continue;
    end
    S = load(fn);
    % Basic sanity: ensure ChannelFlag matches F rows
    if isfield(S,'ChannelFlag') && numel(S.ChannelFlag) ~= size(S.F,1)
        % try to pad flags with 1s or trim to match F
        if numel(S.ChannelFlag) < size(S.F,1)
            S.ChannelFlag = [S.ChannelFlag(:); ones(size(S.F,1)-numel(S.ChannelFlag),1)];
        else
            S.ChannelFlag = S.ChannelFlag(1:size(S.F,1));
        end
    end
    Flist{end+1} = S.F; %#ok<AGROW>
    Slast = S;
end

if isempty(Flist)
    outFile = '';
    fprintf('  No trials found to average for Run %d (%s)\n', runNo, tag);
    return;
end

% Align sample counts if needed (trim to min)
nSamp = min(cellfun(@(x) size(x,2), Flist));
if ~all(cellfun(@(x) size(x,2)==nSamp, Flist))
    for k = 1:numel(Flist)
        Flist{k} = Flist{k}(:,1:nSamp);
    end
    if isfield(Slast,'Time') && numel(Slast.Time) >= nSamp
        Slast.Time = Slast.Time(1:nSamp);
    end
end

% Ensure same #channels across files
nChan = min(cellfun(@(x) size(x,1), Flist));
if ~all(cellfun(@(x) size(x,1)==nChan, Flist))
    for k = 1:numel(Flist)
        Flist{k} = Flist{k}(1:nChan, :);
    end
    if isfield(Slast,'ChannelFlag') && numel(Slast.ChannelFlag) >= nChan
        Slast.ChannelFlag = Slast.ChannelFlag(1:nChan);
    end
end

% Compute average
Favg = mean(cat(3, Flist{:}), 3);

Savg = Slast;
Savg.F       = Favg;
Savg.nAvg    = numel(Flist);
Savg.Comment = sprintf('Avg |Run=%d| |%s| (N=%d) [from T]', runNo, tag, Savg.nAvg);
if isfield(Savg,'History') && iscell(Savg.History)
    Savg.History(end+1,1:3) = {datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
        'avg_by_rt_from_table', sprintf('Run=%d %s N=%d', runNo, tag, Savg.nAvg)};
end

outFile = fullfile(runDir, sprintf('data_%d_avg_%s.mat', runNo, tag));
save(outFile, '-struct', 'Savg');
fprintf('RUN %d: wrote %s (trials=%d, missing=%d)\n', runNo, outFile, numel(Flist), missing);
end

function [A,B] = align_FS(A,B)
% Align channels/samples between two avg structs by trimming to common min
nSamp = min(size(A.F,2), size(B.F,2));
nChan = min(size(A.F,1), size(B.F,1));
A.F = A.F(1:nChan,1:nSamp);
B.F = B.F(1:nChan,1:nSamp);
if isfield(A,'Time'), A.Time = A.Time(1:nSamp); end
if isfield(B,'Time'), B.Time = B.Time(1:nSamp); end
if isfield(A,'ChannelFlag'), A.ChannelFlag = A.ChannelFlag(1:nChan); end
if isfield(B,'ChannelFlag'), B.ChannelFlag = B.ChannelFlag(1:nChan); end
end
