function out = collect_trials_to_append(T, subjectDir, varargin)
% Collect Brainstorm trials into subjectDir/Append/<group>/ as standalone trials.
% Assumes per-run trial files named: data_2_trialNNN.mat in each *_Run_XX folder.
%
% REQUIRED table columns (case-insensitive): run, trial, rt_category
% rt_category can be numeric (1/2/3) or strings like "fast"/"slow"/"timeout".

% -------- Options --------
p = inputParser;
addParameter(p,'FilePrefix','data_2_',@(s)ischar(s)||isstring(s));
addParameter(p,'Groups',[ ...
    struct('name','rt1_fast','gid',1,'rtcat',1), ...
    struct('name','rt2_slow','gid',2,'rtcat',2), ...
    struct('name','rt3_timeout','gid',3,'rtcat',3)], @isstruct);
parse(p,varargin{:});
opt = p.Results;

% Normalize Groups to a struct array
opt.Groups = normalize_groups(opt.Groups);
nG = numel(opt.Groups);

% -------- Setup --------
appendDir = fullfile(subjectDir, 'Append');
if ~exist(appendDir,'dir'), mkdir(appendDir); end

% Pick a channel file template (prefer Run_02 if present)
chanTemplate = pick_channel_template(subjectDir);
if isempty(chanTemplate)
    error('No channel_*.mat file found in run folders.');
end
chanBase = ['channel_' chanTemplate.short];  % the filename well place in each group folder

% -------- Table columns --------
assert(istable(T), 'T must be a table.');
vnRun   = findVar(T, {'run','Run'});
vnTrial = findVar(T, {'trial','Trial','trl','TrialNumber','TrialID'});
vnRT    = findVar(T, {'rt_category','rtcat','rt','RT','category'});

runCol   = T.(vnRun);
trialCol = T.(vnTrial);
rtColRaw = T.(vnRT);

% Convert rt_category to numeric 1/2/3
rtNum = to_rt_numeric(rtColRaw);
assert(numel(rtNum) == height(T), 'rt_category length mismatch with table T.');
rtNum = double(rtNum(:));

% -------- Preallocate output (use cell arrays, not [] ) --------
out = struct( ...
    'group',   {opt.Groups.name}, ...
    'gid',     num2cell([opt.Groups.gid]), ...
    'destDir', repmat({''}, 1, nG), ...
    'nCopied', num2cell(zeros(1,nG)), ...
    'missing', num2cell(zeros(1,nG)), ...
    'files',   repmat({{}}, 1, nG) );

% -------- Process each group --------
for g = 1:nG
    G = opt.Groups(g);
    gDir = fullfile(appendDir, char(G.name));
    if ~exist(gDir,'dir'), mkdir(gDir); end
    out(g).destDir = gDir;

    % Ensure local ChannelFile in this folder
    dstChan = fullfile(gDir, chanBase);
    if ~isfile(dstChan)
        copyfile(chanTemplate.full, dstChan);
    end

    % Row indices for this group (supports scalar or vector G.rtcat)
    rowIdx = find(ismember(rtNum, G.rtcat));
    nCopied = 0; nMissing = 0; dstList = {};

    for ii = 1:numel(rowIdx)
        i = rowIdx(ii);
        runNo   = numify(runCol(i));
        trialNo = numify(trialCol(i));

        runDir = find_run_dir(subjectDir, runNo);
        if isempty(runDir)
            fprintf('RUN %d: folder not found.\n', runNo);
            nMissing = nMissing + 1;
            continue;
        end

        src = fullfile(runDir, sprintf('%strial%03d.mat', opt.FilePrefix, trialNo));
        if ~isfile(src)
            fprintf('  MISSING: %s\n', src);
            nMissing = nMissing + 1;
            continue;
        end

        % Load and sanitize
        S = load(src);

        % Fix ChannelFlag vs F rows if needed
        if isfield(S,'ChannelFlag') && numel(S.ChannelFlag) ~= size(S.F,1)
            if numel(S.ChannelFlag) < size(S.F,1)
                S.ChannelFlag = [S.ChannelFlag(:); ones(size(S.F,1)-numel(S.ChannelFlag),1)];
            else
                S.ChannelFlag = S.ChannelFlag(1:size(S.F,1));
            end
        end

        % Local channel file reference
        S.ChannelFile = chanBase;

        % Tag comment and history
        tag = sprintf('|Group=%d %s| |Run=%02d| |Trial=%03d|', G.gid, G.name, runNo, trialNo);
        if ~isfield(S,'Comment') || isempty(S.Comment), S.Comment = ''; end
        if ~contains(S.Comment, tag), S.Comment = strtrim([S.Comment ' ' tag]); end
        if isfield(S,'History') && iscell(S.History)
            S.History(end+1,1:3) = {datestr(now,'yyyy-mm-dd HH:MM:SS'), 'collect_trials_to_append', tag};
        end

        % Save to destination
        dst = fullfile(gDir, sprintf('%sgrp%d_run%02d_trial%03d.mat', opt.FilePrefix, G.gid, runNo, trialNo));
        save(dst, '-struct', 'S');
        nCopied = nCopied + 1;
        dstList{end+1} = dst; %#ok<AGROW>
    end

    out(g).nCopied = nCopied;
    out(g).missing = nMissing;
    out(g).files   = dstList;
    fprintf('Group %s: copied %d, missing %d -> %s\n', G.name, nCopied, nMissing, gDir);
end

fprintf('\nDone. In Brainstorm: right-click subject ? Reload to see Append/* folders.\n');

end

% ================= Helpers =================

function groups = normalize_groups(G)
% Accept scalar struct with vector fields OR proper struct array.
if numel(G) > 1, groups = G; return; end
n = max([numel(G), numel(G.rtcat), numel(G.gid), numel(G.name)]);
if n == 1, groups = G; return; end
names = G.name; if ~iscell(names), names = cellstr(string(names)); end
gids  = G.gid(:);   if numel(gids) < n, gids(end+1:n) = gids(end); end
rtcs  = G.rtcat(:); if numel(rtcs) < n, rtcs(end+1:n) = rtcs(end); end
groups = repmat(struct('name','','gid',[],'rtcat',[]), 1, n);
for k = 1:n
    groups(k).name  = string(names{k});
    groups(k).gid   = gids(k);
    groups(k).rtcat = rtcs(k);
end
end

function runDir = find_run_dir(subjectDir, runNo)
d = dir(fullfile(subjectDir, sprintf('*_Run_%02d', runNo)));
if isempty(d), runDir = ''; else, runDir = fullfile(d(1).folder, d(1).name); end
end

function out = pick_channel_template(subjectDir)
% Prefer Run_02; else first run containing a channel_*.mat
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
rtNum = str2double(s); % "1","2","3"
L = lower(strtrim(s));
rtNum(ismissing(rtNum) & L=="fast")              = 1;
rtNum(ismissing(rtNum) & L=="slow")              = 2;
rtNum(ismissing(rtNum) & contains(L,'timeout'))  = 3;
rtNum(isnan(rtNum)) = 0; % non-matching labels drop out
end
