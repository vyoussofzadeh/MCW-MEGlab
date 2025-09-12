function bst_bucket_trials_by_response(T_or_path, subjectName, runOffset, varargin)
% BST_BUCKET_TRIALS_BY_RESPONSE
% Move trial epochs from per-run *_resp conditions into aggregate Resp1/Resp2
% conditions using trial-response assignments provided in a table T.
%
% INPUTS
%   T_or_path   : table or path to CSV with columns: run, trial_idx, resp
%                 - run: string like "cidmeg004_times_run1"
%                 - trial_idx: original trial index within the run (1..N)
%                 - resp: one of {"resp1","resp2","timeout"} (case-insensitive)
%   subjectName : Brainstorm subject name (string/char)
%   runOffset   : integer offset to add to parsed run numbers (e.g., 1 if Run_01 is empty-room)
%
% NAME-VALUE OPTIONS
%   'SourceSuffix' : char/string, default '_resp'
%   'TargetResp1'  : char/string, default 'Resp1'
%   'TargetResp2'  : char/string, default 'Resp2'
%   'DryRun'       : logical,     default false
%   'DoChecks'     : logical,     default true  (sanity-check T before moving)
%   'MapMode'      : 'auto'|'rank'|'number' (default 'auto')
%       - 'rank'   : pair by rank (ascending order) of T.trial_idx (resp-only) to the file list.
%                    This is robust when Brainstorm renumbers epochs sequentially (1..nResp).
%       - 'number' : match T.trial_idx to the numeric suffix in filenames '...trialNNN.mat'.
%                    Use only if filenames preserve original trial numbers.
%       - 'auto'   : use 'rank' if filenames look sequential (1..N), otherwise try 'number'.
%   'PreviewN'    : integer,      default 10     (how many example moves to print in DryRun)
%   'Verbose'     : logical,      default true
%
% EXAMPLE
%   bst_bucket_trials_by_response('cidmeg004_times.csv','cidmeg004',0,'DryRun',true);
%
% REQUIREMENTS
%   - Brainstorm GUI/session with database loaded (uses bst_get, bst_process)

% -----------------------------
% Parse inputs
% -----------------------------
p = inputParser;
addRequired(p,'T_or_path');
addRequired(p,'subjectName',@(s)ischar(s)||isstring(s));
addRequired(p,'runOffset',@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'SourceSuffix','_resp',@(s)ischar(s)||isstring(s));
addParameter(p,'TargetResp1','Resp1',@(s)ischar(s)||isstring(s));
addParameter(p,'TargetResp2','Resp2',@(s)ischar(s)||isstring(s));
addParameter(p,'DryRun',false,@islogical);
addParameter(p,'DoChecks',true,@islogical);
addParameter(p,'MapMode','auto',@(s)any(strcmpi(s,{'auto','rank','number'})));
addParameter(p,'PreviewN',10,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'Verbose',true,@islogical);
parse(p,T_or_path,subjectName,runOffset,varargin{:});
opt = p.Results;
opt.MapMode = lower(string(opt.MapMode));

% -----------------------------
% Load / normalize T
% -----------------------------
if istable(opt.T_or_path)
    T = opt.T_or_path;
else
    T = readtable(char(opt.T_or_path));
end
reqCols = {'run','trial_idx','resp'};
assert(all(ismember(reqCols, T.Properties.VariableNames)), ...
    'T must have columns: run, trial_idx, resp');
T.run  = string(T.run);
T.resp = lower(string(T.resp));

% -----------------------------
% Brainstorm subject exists?
% -----------------------------
[sSubject, iSubject] = bst_get('Subject', char(opt.subjectName), 1); %#ok<ASGLU>
assert(~isempty(iSubject), 'Subject not found in Brainstorm: %s', char(opt.subjectName));

% -----------------------------
% Optional sanity checks
% -----------------------------
if opt.DoChecks
    run_sanity_checks(T);
end

% -----------------------------
% Iterate runs in T (stable)
% -----------------------------
runs_csv = unique(T.run, 'stable');
toMove_R1 = {}; toMove_R2 = {};
movedInfoR1 = {}; movedInfoR2 = {};  % for DryRun preview
skipped = 0;

for k = 1:numel(runs_csv)
    runCsv = runs_csv(k);
    rnum   = local_parse_trailing_int(runCsv);     % "...run9" -> 9
    bsRun  = rnum + opt.runOffset;                 % apply offset
    condName   = sprintf('%s_Run_%02d%s', char(opt.subjectName), bsRun, char(opt.SourceSuffix));
    fullCondId = sprintf('%s/%s', char(opt.subjectName), condName);

    % Fetch condition
    [~, iStudy] = bst_get('StudyWithCondition', fullCondId);
    if isempty(iStudy)
        warning('Condition not found: %s', fullCondId);
        continue;
    end
    sStudy = bst_get('Study', iStudy);

    % Collect data file basenames
    files = {};
    if isfield(sStudy,'Data') && ~isempty(sStudy.Data)
        files = {sStudy.Data.FileName};
    end
    if isempty(files)
        warning('No data files in condition: %s', condName);
        continue;
    end
    % Basenames
    [~, names, exts] = cellfun(@fileparts, files, 'UniformOutput', false);
    bases = strcat(names, exts);

    % Keep only files ending with trialNNN.mat (generic)
    isData = ~cellfun('isempty', regexp(bases, 'trial(\d+)\.mat$', 'once'));
    files  = files(isData);
    bases  = bases(isData);

    if isempty(files)
        warning('No files matching pattern "*trialNNN.mat" in: %s', condName);
        continue;
    end

    % Extract trial numbers (if present)
    trialNums = nan(numel(bases),1);
    for j = 1:numel(bases)
        tok = regexp(bases{j}, 'trial(\d+)\.mat$', 'tokens', 'once');
        if ~isempty(tok)
            trialNums(j) = str2double(tok{1});
        end
    end

    % Sort by trialNums when available; else keep original order
    haveNums = ~isnan(trialNums);
    if any(haveNums)
        [trialNums, ord] = sort(trialNums(haveNums));
        filesSorted      = files(haveNums);
        filesSorted      = filesSorted(ord);
        trialNumsSorted  = trialNums(:);
    else
        filesSorted     = files;    % fallback
        trialNumsSorted = (1:numel(filesSorted)).'; %#ok<NASGU>  % placeholder
    end

    % Subset T rows for this run, RESP-ONLY
    selRun  = T.run == runCsv & ismember(T.resp, ["resp1","resp2"]);
    tIdx    = T.trial_idx(selRun);
    tRsp    = T.resp(selRun);

    if isempty(tIdx)
        if opt.Verbose
            fprintf('%s: No resp1/resp2 rows in T for this run (skipping).\n', condName);
        end
        continue;
    end

    % Decide mapping mode for this run
    mapModeRun = opt.MapMode;
    if mapModeRun == "auto"
        % If filenames look strictly sequential 1..N, prefer rank mapping.
        if any(haveNums) && isequal((1:numel(filesSorted)).', trialNumsSorted)
            mapModeRun = "rank";
        else
            % Otherwise try number mapping if any overlap; else use rank.
            if any(haveNums) && any(ismember(tIdx(:), trialNumsSorted))
                mapModeRun = "number";
            else
                mapModeRun = "rank";
            end
        end
    end

    switch mapModeRun
        case "rank"
            % Pair by rank: sort T indices asc; assign j-th responded trial to j-th file.
            [tIdxSorted, ordT] = sort(tIdx(:)); %#ok<NASGU>
            tRspSorted         = tRsp(ordT);
            nFiles = numel(filesSorted);
            nRows  = numel(tRspSorted);
            nPair  = min(nFiles, nRows);

            if nFiles ~= nRows
                warning('%s: %d files vs %d T rows; pairing by rank (min=%d).', ...
                        condName, nFiles, nRows, nPair);
            end

            for j = 1:nPair
                f = filesSorted{j};
                if tRspSorted(j) == "resp1"
                    toMove_R1{end+1} = f; %#ok<AGROW>
                    if opt.DryRun, movedInfoR1{end+1} = sprintf('%s: #%d -> Resp1', condName, j); end %#ok<AGROW>
                else
                    toMove_R2{end+1} = f; %#ok<AGROW>
                    if opt.DryRun, movedInfoR2{end+1} = sprintf('%s: #%d -> Resp2', condName, j); end %#ok<AGROW>
                end
            end
            skipped = skipped + abs(nFiles - nRows);

        case "number"
            % Map by original trial number found in filename suffix trialNNN
            if ~any(haveNums)
                warning('%s: No numeric trial suffixes found; falling back to rank mapping.', condName);
                % Fallback to rank
                [tIdxSorted, ordT] = sort(tIdx(:)); %#ok<NASGU>
                tRspSorted         = tRsp(ordT);
                nFiles = numel(filesSorted);
                nRows  = numel(tRspSorted);
                nPair  = min(nFiles, nRows);
                for j = 1:nPair
                    f = filesSorted{j};
                    if tRspSorted(j) == "resp1"
                        toMove_R1{end+1} = f; %#ok<AGROW>
                        if opt.DryRun, movedInfoR1{end+1} = sprintf('%s: #%d -> Resp1', condName, j); end %#ok<AGROW>
                    else
                        toMove_R2{end+1} = f; %#ok<AGROW>
                        if opt.DryRun, movedInfoR2{end+1} = sprintf('%s: #%d -> Resp2', condName, j); end %#ok<AGROW>
                    end
                end
                skipped = skipped + abs(nFiles - nRows);
            else
                % Build finder
                trialNumsRun = trialNumsSorted; % already sorted with filesSorted
                for j = 1:numel(tIdx)
                    pos = find(trialNumsRun == tIdx(j), 1, 'first');
                    if isempty(pos)
                        skipped = skipped + 1;
                        continue;
                    end
                    f = filesSorted{pos};
                    if tRsp(j) == "resp1"
                        toMove_R1{end+1} = f; %#ok<AGROW>
                        if opt.DryRun, movedInfoR1{end+1} = sprintf('%s: trial%03d -> Resp1', condName, tIdx(j)); end %#ok<AGROW>
                    else
                        toMove_R2{end+1} = f; %#ok<AGROW>
                        if opt.DryRun, movedInfoR2{end+1} = sprintf('%s: trial%03d -> Resp2', condName, tIdx(j)); end %#ok<AGROW>
                    end
                end
            end

        otherwise
            error('Unknown MapMode: %s', mapModeRun);
    end
end

% -----------------------------
% Dry run preview or execute
% -----------------------------
if opt.DryRun
    fprintf('[DryRun] Would move %d files to "%s" and %d files to "%s". Skipped=%d.\n', ...
        numel(toMove_R1), char(opt.TargetResp1), numel(toMove_R2), char(opt.TargetResp2), skipped);

    if opt.Verbose && opt.PreviewN > 0
        n1 = min(opt.PreviewN, numel(movedInfoR1));
        n2 = min(opt.PreviewN, numel(movedInfoR2));
        if n1>0
            fprintf('  Examples -> %s:\n', char(opt.TargetResp1));
            fprintf('    %s\n', movedInfoR1{1:n1});
        end
        if n2>0
            fprintf('  Examples -> %s:\n', char(opt.TargetResp2));
            fprintf('    %s\n', movedInfoR2{1:n2});
        end
    end
    return;
end

movedR1 = 0; movedR2 = 0;
if ~isempty(toMove_R1)
    bst_process('CallProcess', 'process_set_condition', toMove_R1, [], ...
        'subjectname', char(opt.subjectName), ...
        'condition',   char(opt.TargetResp1), ...
        'isabsolute',  0, 'createcond', 1);
    movedR1 = numel(toMove_R1);
end
if ~isempty(toMove_R2)
    bst_process('CallProcess', 'process_set_condition', toMove_R2, [], ...
        'subjectname', char(opt.subjectName), ...
        'condition',   char(opt.TargetResp2), ...
        'isabsolute',  0, 'createcond', 1);
    movedR2 = numel(toMove_R2);
end

fprintf('Moved %d files to "%s" and %d files to "%s". Skipped=%d.\n', ...
    movedR1, char(opt.TargetResp1), movedR2, char(opt.TargetResp2), skipped);

end % main function

% =============================
% Helpers
% =============================
function n = local_parse_trailing_int(s)
m = regexp(char(s), '(\d+)$', 'tokens', 'once');
assert(~isempty(m), 'Cannot parse run number from "%s"', s);
n = str2double(m{1});
end

function run_sanity_checks(T)
% Basic validation of T prior to moving files.
assert(all(ismember(T.resp, ["resp1","resp2","timeout"])), ...
    'Unexpected values in T.resp. Expected {"resp1","resp2","timeout"}');

% No duplicate (run, trial_idx) among resp1/resp2
T2 = T(ismember(T.resp,["resp1","resp2"]), {'run','trial_idx'});
if ~isempty(T2)
    % Use table grouping if available
    try
        G = groupsummary(T2, {'run','trial_idx'}, 'numel');
        assert(all(G.GroupCount==1), 'Duplicate (run,trial_idx) found among resp1/resp2 rows in T.');
    catch
        % Fallback without groupsummary
        key = strcat(T2.run, "__", string(T2.trial_idx));
        [u,~,g] = unique(key, 'stable');
        counts = accumarray(g, 1);
        if any(counts>1)
            bad = u(counts>1);
            error('Duplicate (run,trial_idx) rows in T for: %s', strjoin(bad(1:min(5,end)), ', '));
        end
    end
end
end
