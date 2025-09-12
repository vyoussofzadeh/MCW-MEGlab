function out = collect_trials_to_groups_acc(T, subjectDir, varargin)
% COLLECT_TRIALS_TO_GROUPS_ACC
% Split trials into FAST/SLOW/TIMEOUT and further into CORRECT/INCORRECT.
% Saves the underlying trial .mat files into group folders:
%   fast_correct, fast_incorrect, slow_correct, slow_incorrect, timeout
%
% out = collect_trials_to_groups_acc(T, subjectDir, 'FilePrefix','data_2_');
%
% REQUIRED:
%   T           table with trial metadata. Expected columns (default names):
%               - GroupVar: 'group' or 'resp_group' or numeric (1=slow, 2=fast, 3=timeout)
%               - AccVar:   logical/0-1 or text ('correct'/'incorrect')
%               - FileVar:  path/filename of each trial file (if absent, we build from Run/Trial)
%               Optional: 'Run','Trial' to build filenames if FileVar missing
%   subjectDir  Brainstorm "data/<subject>" directory that contains the trial .mat files
%
% OPTIONS (Name-Value):
%   'FilePrefix'   (char) default 'data_2_'
%   'GroupVar'     (char) default auto-detect: first of {'group','resp_group','Group','RespGroup'}
%   'AccVar'       (char) default auto-detect: first of {'accuracy','acc','is_correct','Correct'}
%   'FileVar'      (char) default auto-detect: first of {'file','filename','File','FileName','matfile'}
%   'RunVar'       (char) default 'Run'
%   'TrialVar'     (char) default 'Trial'
%   'Action'       (char) 'copy'|'link'|'none'  default 'copy'
%   'MakeDirs'     (logical) create destination folders if missing (default true)
%
% RETURNS:
%   out.groups.<groupname> with fields: dir, n, trials (cellstr)
%
% NOTES:
%   - If Group is numeric, mapping assumed: 1=slow, 2=fast, 3=timeout.
%   - TIMEOUT trials are not split into correct/incorrect.
%
% VY/MCW utility © 2025

% -------- Parse inputs
p = inputParser;
p.addRequired('T', @(x) istable(x));
p.addRequired('subjectDir', @(x) ischar(x) || isstring(x));
p.addParameter('FilePrefix','data_2_', @(x) ischar(x) || isstring(x));
p.addParameter('GroupVar','', @(x) ischar(x) || isstring(x));
p.addParameter('AccVar','', @(x) ischar(x) || isstring(x));
p.addParameter('FileVar','', @(x) ischar(x) || isstring(x));
p.addParameter('RunVar','Run', @(x) ischar(x) || isstring(x));
p.addParameter('TrialVar','Trial', @(x) ischar(x) || isstring(x));
p.addParameter('Action','copy', @(x) any(strcmpi(x,{'copy','link','none'})));
p.addParameter('MakeDirs',true, @islogical);
p.parse(T, subjectDir, varargin{:});
opt = p.Results;

subjectDir = char(subjectDir);
filePrefix = char(opt.FilePrefix);

% -------- Detect column names if not provided
if isempty(opt.GroupVar)
    opt.GroupVar = firstPresentVar(T, {'group','resp_group','Group','RespGroup'});
end
if isempty(opt.AccVar)
    opt.AccVar = firstPresentVar(T, {'accuracy','acc','is_correct','Correct'});
end
if isempty(opt.FileVar)
    opt.FileVar = firstPresentVar(T, {'file','filename','File','FileName','matfile'});
end
assert(~isempty(opt.GroupVar) && any(strcmp(opt.GroupVar, T.Properties.VariableNames)), ...
    'Group column not found in T.');
assert(~isempty(opt.AccVar) && any(strcmp(opt.AccVar, T.Properties.VariableNames)), ...
    'Accuracy column not found in T.');

% -------- Resolve file names
trialFiles = resolveTrialFiles(T, subjectDir, filePrefix, opt.FileVar, char(opt.RunVar), char(opt.TrialVar));

% -------- Build group keys per trial
grpRaw = T.(opt.GroupVar);
accRaw = T.(opt.AccVar);

grpTxt = normalizeGroup(grpRaw);   % 'fast','slow','timeout'
accTF  = normalizeAcc(accRaw);     % true=correct, false=incorrect, NaN=unknown

n = height(T);
keys = strings(n,1);
for i = 1:n
    if strcmp(grpTxt{i}, 'timeout')
        keys(i) = "timeout";
    else
        if ismissing(accTF(i)) || isnan(accTF(i))
            % default unknown accuracy to incorrect for safety, or make separate bucket
            keys(i) = grpTxt{i} + "_incorrect";
        else
            keys(i) = grpTxt{i} + "_" + (accTF(i) * "correct" + (~accTF(i)) * "incorrect");
        end
    end
end

% -------- Buckets
buckets = ["fast_correct","fast_incorrect","slow_correct","slow_incorrect","timeout"];
out = struct; out.groups = struct;

for b = buckets
    idx = strcmp(keys, b);
    if ~any(idx), continue; end

    dstDir = fullfile(subjectDir, char(b));
    if opt.MakeDirs && ~exist(dstDir, 'dir')
        mkdir(dstDir);
    end

    filesB = trialFiles(idx);
    switch lower(opt.Action)
        case 'copy'
            for k = 1:numel(filesB)
                if exist(filesB{k}, 'file')
                    copyfile(filesB{k}, fullfile(dstDir, fileparts_only(filesB{k}, 'name')));
                end
            end
        case 'link'
            % symlink (Unix); on Windows this may require admin. Fallback to copy if fails.
            for k = 1:numel(filesB)
                src = filesB{k};
                dst = fullfile(dstDir, fileparts_only(src, 'name'));
                if exist(src,'file') && ~exist(dst,'file')
                    try
                        system(sprintf('ln -s "%s" "%s"', src, dst));
                    catch
                        copyfile(src, dst);
                    end
                end
            end
        case 'none'
            % do nothing
    end

    % Store in output (valid field names)
    gname = matlab.lang.makeValidName(char(b));
    out.groups.(gname).name   = char(b);
    out.groups.(gname).dir    = dstDir;
    out.groups.(gname).n      = sum(idx);
    out.groups.(gname).trials = filesB;
end

% also return a flat summary
out.summary = struct();
for b = buckets
    f = matlab.lang.makeValidName(char(b));
    if isfield(out.groups, f)
        out.summary.(f) = out.groups.(f).n;
    else
        out.summary.(f) = 0;
    end
end

end % function main

% ================= helpers =================

function v = firstPresentVar(T, candidates)
v = '';
for i = 1:numel(candidates)
    if any(strcmp(candidates{i}, T.Properties.VariableNames))
        v = candidates{i};
        return;
    end
end
end

function files = resolveTrialFiles(T, subjectDir, filePrefix, fileVar, runVar, trialVar)
if ~isempty(fileVar) && any(strcmp(fileVar, T.Properties.VariableNames))
    val = T.(fileVar);
    if iscell(val)
        files = val;
    elseif isstring(val) || ischar(val)
        files = cellstr(val);
    else
        error('Unsupported FileVar type.');
    end
    % make absolute if needed
    for i = 1:numel(files)
        if ~contains(files{i}, filesep)
            files{i} = fullfile(subjectDir, files{i});
        end
    end
else
    % Build from Run/Trial: data_2_run%02d_trial%03d.mat
    assert(all(ismember({runVar,trialVar}, T.Properties.VariableNames)), ...
        'Run/Trial columns not found; cannot build filenames.');
    run  = T.(runVar);
    trl  = T.(trialVar);
    files = cell(height(T),1);
    for i = 1:height(T)
        files{i} = fullfile(subjectDir, sprintf('%srun%02d_trial%03d.mat', filePrefix, run(i), trl(i)));
    end
end
end

function g = normalizeGroup(graw)
% Map numeric -> text; lower-case text; accept synonyms
if isnumeric(graw)
    g = strings(numel(graw),1);
    for i = 1:numel(graw)
        switch graw(i)
            case 1, g(i) = "slow";
            case 2, g(i) = "fast";
            case 3, g(i) = "timeout";
            otherwise, g(i) = "unknown";
        end
    end
elseif iscell(graw) || isstring(graw)
    g = string(graw);
    g = lower(strtrim(g));
    g = replace(g, {'1_slow','2_fast','3_timeout'}, {'slow','fast','timeout'});
else
    error('Unsupported GroupVar type.');
end
end

function acc = normalizeAcc(a)
% Returns true for correct, false for incorrect, NaN for unknown
if islogical(a)
    acc = a;
elseif isnumeric(a)
    acc = a ~= 0;
elseif iscell(a) || isstring(a)
    s = lower(strtrim(string(a)));
    acc = nan(numel(s),1);
    acc(ismember(s, ["correct","hit","true","t","yes","y","1"]))  = true;
    acc(ismember(s, ["incorrect","miss","false","f","no","n","0"])) = false;
else
    acc = nan(size(a));
end
acc = reshape(acc, [], 1);
end

function name = fileparts_only(p, what)
% Return just 'name.ext' of a path
[~, n, e] = fileparts(p);
switch lower(what)
    case 'name', name = [n e];
    otherwise, name = [n e];
end
end