function out = avg_allruns_by_rt_from_table(T, subjectDir, varargin)
% Average Brainstorm trials across ALL runs for one subject using table T.
% T needs columns: run, trial, rt_category (1=fast, 2=slow, 3=timeout).
% Assumes trial files in: subjectDir/*_Run_%02d/data_2_trialNNN.mat
%
% Usage:
% out = avg_allruns_by_rt_from_table(T, '/path/to/mcwa086_v1', ...
%        'IncludeTimeout', false, 'DoContrast', true);

% ---- options ----
p = inputParser;
addParameter(p,'IncludeTimeout',false,@islogical);
addParameter(p,'DoContrast',true,@islogical);
addParameter(p,'Runs',[],@(x)isnumeric(x) || isempty(x));
addParameter(p,'FilePrefix','data_2_',@(s)ischar(s) || isstring(s));  % you said files always start with data_2_
parse(p,varargin{:});
opt = p.Results;

FAST = 1; SLOW = 2; TIMEOUT = 3;

if isempty(opt.Runs)
    runs = unique(T.run(:))';
else
    runs = opt.Runs(:)';
end

runs = runs+1;

% pick save dir (@intra if exists, else first run dir, else subjectDir)
intraDir = fullfile(subjectDir, '@intra');
if ~exist(intraDir, 'dir')
    dFirst = dir(fullfile(subjectDir, sprintf('*_Run_%02d', runs(1))));
    if ~isempty(dFirst), intraDir = fullfile(dFirst(1).folder, dFirst(1).name);
    else, intraDir = subjectDir;
    end
end
if ~exist(intraDir,'dir'), mkdir(intraDir); end

% ----- Gather file lists per category across runs -----
filesFast = {}; filesSlow = {}; filesTO = {};
missFast=0; missSlow=0; missTO=0;

for r = runs
    runDir = find_run_dir(subjectDir, r);
    if isempty(runDir)
        fprintf('RUN %d: folder not found.\n', r);
        continue;
    end
    rows = T(T.run == r, :);
    % fast
    ftr = rows.trial(rows.rt_category == FAST);
    [list, m] = trial_files(runDir, opt.FilePrefix, ftr);
    filesFast = [filesFast list]; %#ok<AGROW>
    missFast = missFast + m;
    % slow
    str = rows.trial(rows.rt_category == SLOW);
    [list, m] = trial_files(runDir, opt.FilePrefix, str);
    filesSlow = [filesSlow list]; %#ok<AGROW>
    missSlow = missSlow + m;
    % timeout
    if opt.IncludeTimeout
        ttr = rows.trial(rows.rt_category == TIMEOUT);
        [list, m] = trial_files(runDir, opt.FilePrefix, ttr);
        filesTO = [filesTO list]; %#ok<AGROW>
        missTO = missTO + m;
    end
end

% ----- Average & save -----
out = struct('avgFastFile','', 'avgSlowFile','', 'avgTOFile','', 'contrastFile','', ...
             'nFast',numel(filesFast), 'nSlow',numel(filesSlow), 'nTimeout',numel(filesTO), ...
             'missFast',missFast, 'missSlow',missSlow, 'missTO',missTO);

if ~isempty(filesFast)
    out.avgFastFile = do_avg(filesFast, intraDir, sprintf('%savg_AllRuns_RTfast_fromT.mat', opt.FilePrefix));
end
if ~isempty(filesSlow)
    out.avgSlowFile = do_avg(filesSlow, intraDir, sprintf('%savg_AllRuns_RTslow_fromT.mat', opt.FilePrefix));
end
if ~isempty(filesTO)
    out.avgTOFile = do_avg(filesTO, intraDir, sprintf('%savg_AllRuns_RTtimeout_fromT.mat', opt.FilePrefix));
end

% optional contrast: Fast - Slow
if opt.DoContrast && ~isempty(out.avgFastFile) && ~isempty(out.avgSlowFile)
    Sfast = load(out.avgFastFile);
    Sslow = load(out.avgSlowFile);
    [Sfast, Sslow] = align_FS(Sfast, Sslow);
    Sdiff = Sfast;
    Sdiff.F = Sfast.F - Sslow.F;
    Sdiff.Comment = 'Contrast |AllRuns| |RT=fast-minus-slow| [from T]';
    if isfield(Sdiff,'History') && iscell(Sdiff.History)
        Sdiff.History(end+1,1:3) = {datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
                                    'avg_allruns_by_rt_from_table','Fast - Slow (AllRuns)'};
    end
    out.contrastFile = fullfile(intraDir, sprintf('%scontrast_AllRuns_RTfast_minus_RTslow_fromT.mat', opt.FilePrefix));
    save(out.contrastFile, '-struct', 'Sdiff');
    fprintf('Wrote %s\n', out.contrastFile);
end

% Summary
fprintf('\nSUMMARY (AllRuns): Fast=%d (miss %d), Slow=%d (miss %d), Timeout=%d (miss %d)\n', ...
    out.nFast, out.missFast, out.nSlow, out.missSlow, out.nTimeout, out.missTO);
fprintf('Saved to: %s\n', intraDir);

end

% --------- helpers ---------
function runDir = find_run_dir(subjectDir, runNo)
pat = sprintf('*_Run_%02d', runNo);
d = dir(fullfile(subjectDir, pat));
if isempty(d), runDir = ''; else, runDir = fullfile(d(1).folder, d(1).name); end
end

function [flist, nMissing] = trial_files(runDir, prefix, trials)
flist = {};
nMissing = 0;
for t = trials(:)'
    fn = fullfile(runDir, sprintf('%strial%03d.mat', prefix, t));
    if isfile(fn), flist{end+1} = fn; else, nMissing = nMissing + 1; fprintf('  MISSING: %s\n', fn); end %#ok<AGROW>
end
end

function outFile = do_avg(fileList, saveDir, outName)
% Load all, fix flag sizes, align dims, average, write MAT
Flist = {};
lastS = [];
for i = 1:numel(fileList)
    S = load(fileList{i});
    % Fix ChannelFlag length vs F rows if needed
    if isfield(S,'ChannelFlag') && numel(S.ChannelFlag) ~= size(S.F,1)
        if numel(S.ChannelFlag) < size(S.F,1)
            S.ChannelFlag = [S.ChannelFlag(:); ones(size(S.F,1)-numel(S.ChannelFlag),1)];
        else
            S.ChannelFlag = S.ChannelFlag(1:size(S.F,1));
        end
    end
    Flist{end+1} = S.F(1:306,:); %#ok<AGROW>
    lastS = S;
end
% Align samples/channels to min common
nSamp = min(cellfun(@(x) size(x,2), Flist));
nChan = min(cellfun(@(x) size(x,1), Flist));
for k = 1:numel(Flist), Flist{k} = Flist{k}(1:nChan,1:nSamp); end
if isfield(lastS,'Time') && numel(lastS.Time) >= nSamp, lastS.Time = lastS.Time(1:nSamp); end
if isfield(lastS,'ChannelFlag') && numel(lastS.ChannelFlag) >= nChan, lastS.ChannelFlag = lastS.ChannelFlag(1:nChan); end

Favg = mean(cat(3, Flist{:}), 3);
Savg = lastS;
Savg.F = Favg;
Savg.nAvg = numel(Flist);
Savg.Comment = sprintf('Avg |AllRuns| (N=%d) [from T]', Savg.nAvg);
if isfield(Savg,'History') && iscell(Savg.History)
    Savg.History(end+1,1:3) = {datestr(now,'yyyy-mm-dd HH:MM:SS'), 'avg_allruns_by_rt_from_table', outName};
end

outFile = fullfile(saveDir, outName);
save(outFile, '-struct', 'Savg');
fprintf('Wrote %s (trials=%d)\n', outFile, numel(Flist));
end

function [A,B] = align_FS(A,B)
nSamp = min(size(A.F,2), size(B.F,2));
nChan = min(size(A.F,1), size(B.F,1));
A.F = A.F(1:nChan,1:nSamp);
B.F = B.F(1:nChan,1:nSamp);
if isfield(A,'Time'), A.Time = A.Time(1:nSamp); end
if isfield(B,'Time'), B.Time = B.Time(1:nSamp); end
if isfield(A,'ChannelFlag'), A.ChannelFlag = A.ChannelFlag(1:nChan); end
if isfield(B,'ChannelFlag'), B.ChannelFlag = B.ChannelFlag(1:nChan); end
end
