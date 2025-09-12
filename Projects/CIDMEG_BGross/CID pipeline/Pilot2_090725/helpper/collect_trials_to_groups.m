function out = collect_trials_to_groups(T, subjectDir, varargin)
% Copy Brainstorm trials into subjectDir/<group>/ as standalone trials.
% REQUIRED table columns (case-insensitive): run, trial, rt_category
% rt_category may be numeric 1/2/3 or strings like "fast"/"slow"/"timeout".
%
% Default groups (per your request):
%  - 1_slow     <- rt_category == 2
%  - 2_fast     <- rt_category == 1
%  - 3_timeout  <- rt_category == 3
%
% Example:
% out = collect_trials_to_groups(T, '/.../mcwa086_v1', 'FilePrefix','data_2_');

% -------- Options --------
p = inputParser;
addParameter(p,'FilePrefix','data_2_',@(s)ischar(s)||isstring(s));
% Folder names and mapping (note: labels vs categories are intentionally swapped)
addParameter(p,'Groups',[ ...
    struct('name','slow',    'gid',1,'rtcat',2), ...
    struct('name','fast',    'gid',2,'rtcat',1), ...
    struct('name','timeout', 'gid',3,'rtcat',3)], @isstruct);
parse(p,varargin{:});
opt = p.Results;

% Normalize Groups to a clean 1xN struct array
opt.Groups = normalize_groups(opt.Groups);
nG = numel(opt.Groups);

% -------- Channel template (prefer Run_02) --------
chanTemplate = pick_channel_template(subjectDir);
if isempty(chanTemplate)
    error('No channel_*.mat found under %s/*_Run_*', subjectDir);
end
chanBase = ['channel_' chanTemplate.short];  % we will place this filename in each destination folder

% -------- Table columns --------
assert(istable(T), 'T must be a table.');
vnRun   = findVar(T, {'run','Run'});
vnTrial = findVar(T, {'trial','Trial','trl','TrialNumber','TrialID'});
vnRT    = findVar(T, {'rt_category','rtcat','rt','RT','category'});

runCol   = T.(vnRun);
trialCol = T.(vnTrial);
rtColRaw = T.(vnRT);

% Make rt_category numeric 1/2/3
rtNum = to_rt_numeric(rtColRaw);
assert(numel(rtNum) == height(T), 'rt_category length mismatch with T.');
rtNum = double(rtNum(:));

% -------- Preallocate output (column shapes) --------
groupNames = reshape(cellstr(string({opt.Groups.name})), [], 1);
gids       = num2cell(reshape([opt.Groups.gid], [], 1));
destDirs   = repmat({''}, nG, 1);
nCopiedC   = num2cell(zeros(nG,1));
missingC   = num2cell(zeros(nG,1));
filesC     = repmat({{}}, nG, 1);

out = struct('group',groupNames,'gid',gids,'destDir',destDirs, ...
             'nCopied',nCopiedC,'missing',missingC,'files',filesC);

% -------- Copy per group --------
for g = 1:nG
    G = opt.Groups(g);
    gDir = fullfile(subjectDir, char(G.name));
    if ~exist(gDir,'dir'), mkdir(gDir); end
    out(g).destDir = gDir;

    % Ensure a local ChannelFile copy
    dstChan = fullfile(gDir, chanBase);
    if ~isfile(dstChan)
        copyfile(chanTemplate.full, dstChan);
    end

    % Which rows belong to this group?
    rowIdx = find(ismember(rtNum, G.rtcat));   % supports scalar or vector G.rtcat
    nCopied = 0; nMissing = 0; dstList = {};

    for ii = 1:numel(rowIdx)
        i = rowIdx(ii);
        runNo   = numify(runCol(i));
        trialNo = numify(trialCol(i));

        runDir = find_run_dir(subjectDir, runNo);
        if isempty(runDir)
            fprintf('RUN %d: folder not found.\n', runNo);
            nMissing = nMissing + 1; continue;
        end

        src = fullfile(runDir, sprintf('%strial%03d.mat', opt.FilePrefix, trialNo));
        if ~isfile(src)
            fprintf('  MISSING: %s\n', src);
            nMissing = nMissing + 1; continue;
        end

        % Load and sanitize the trial
        S = load(src);

        % Fix ChannelFlag size if needed
        if isfield(S,'ChannelFlag') && numel(S.ChannelFlag) ~= size(S.F,1)
            if numel(S.ChannelFlag) < size(S.F,1)
                S.ChannelFlag = [S.ChannelFlag(:); ones(size(S.F,1)-numel(S.ChannelFlag),1)];
            else
                S.ChannelFlag = S.ChannelFlag(1:size(S.F,1));
            end
        end

        % Point to the local channel file
        S.ChannelFile = chanBase;

        % Tag comment + history
        tag = sprintf('|Group=%d %s| |Run=%02d| |Trial=%03d|', G.gid, G.name, runNo, trialNo);
        if ~isfield(S,'Comment') || isempty(S.Comment), S.Comment = ''; end
        if ~contains(S.Comment, tag), S.Comment = strtrim([S.Comment ' ' tag]); end
        if isfield(S,'History') && iscell(S.History)
            S.History(end+1,1:3) = {datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
                                    'collect_trials_to_groups', tag};
        end

        % Destination filename (unique within group)
        dst = fullfile(gDir, sprintf('%srun%02d_trial%03d.mat', opt.FilePrefix, runNo, trialNo));
        save(dst, '-struct', 'S');
        nCopied = nCopied + 1;
        dstList{end+1} = dst; %#ok<AGROW>
    end

    out(g).nCopied = nCopied;
    out(g).missing = nMissing;
    out(g).files   = dstList;
    fprintf('Group %s: copied %d, missing %d -> %s\n', G.name, nCopied, nMissing, gDir);
end

fprintf('\nDone. In Brainstorm: right-click the subject ? Reload to see 1_slow / 2_fast / 3_timeout.\n');

end

% ================= Helpers =================

function groups = normalize_groups(G)
% Accept scalar struct with vector fields OR proper struct array; flatten if needed.
if isscalar(G)
    n = max([numel(G.name), numel(G.gid), numel(G.rtcat), 1]);
    names = G.name; if ~iscell(names), names = cellstr(string(names)); end
    gids  = G.gid;   if ~iscell(gids),  gids  = num2cell(gids(:).');  end
    rtcs  = G.rtcat; if ~iscell(rtcs),  rtcs  = num2cell(rtcs(:).');  end
    names = names(:).'; if numel(names) < n, names(end+1:n) = names(end); end
    gids  = gids(:).';  if numel(gids)  < n, gids(end+1:n)  = gids(end);  end
    rtcs  = rtcs(:).';  if numel(rtcs)  < n, rtcs(end+1:n)  = rtcs(end);  end
    groups = repmat(struct('name','','gid',[],'rtcat',[]), 1, n);
    for k = 1:n
        groups(k).name  = string(names{k});
        groups(k).gid   = gids{k};
        groups(k).rtcat = rtcs{k};
    end
else
    G = G(:);
    groups = reshape(G, 1, []); % row vector
end
end

function runDir = find_run_dir(subjectDir, runNo)
d = dir(fullfile(subjectDir, sprintf('*_Run_%02d', runNo)));
if isempty(d), runDir = ''; else, runDir = fullfile(d(1).folder, d(1).name); end
end

function out = pick_channel_template(subjectDir)
% Prefer Run_02; else first run with a channel_*.mat
out = struct('full','','short','');
d = dir(fullfile(subjectDir, '*_Run_*'));
if isempty(d), return; end
pref = find(contains({d.name}, '_Run_02'), 1, 'first');
order = [pref, setdiff(1:numel(d), pref)];
order(isnan(order)) = [];
for k = order
    dd = dir(fullfile(d(k).folder, d(k).name, 'channel_*.mat'));
    if ~isempty(dd)
        out.full  = fullfile(dd(1).folder, dd(1).name);
        out.short = erase(dd(1).name, 'channel_');
        return;
    end
end
end

function name = findVar(T, candidates)
names = T.Properties.VariableNames;
for c = 1:numel(candidates)
    ix = find(strcmpi(names, candidates{c}), 1);
    if ~isempty(ix), name = names{ix}; return; end
end
error('Could not find any of these columns in T: %s', strjoin(candidates, ', '));
end

function x = numify(v)
if iscell(v), v = v{1}; end
if iscategorical(v), v = string(v); end
if isstring(v) || ischar(v), x = str2double(string(v));
else, x = double(v);
end
end

function rtNum = to_rt_numeric(rtCol)
if isnumeric(rtCol), rtNum = rtCol; return; end
s = string(rtCol);
rtNum = str2double(s);                 % handles "1","2","3"
L = lower(strtrim(s));
rtNum(ismissing(rtNum) & L=="fast")              = 1;
rtNum(ismissing(rtNum) & L=="slow")              = 2;
rtNum(ismissing(rtNum) & contains(L,'timeout'))  = 3;
rtNum(isnan(rtNum)) = 0; % non-matching labels drop out
end
