function bst_bucket_trials_by_response_copy(T_or_path, subjectName, runOffset, varargin)
% Copy trials into Resp1/Resp2 WITHOUT altering original per-run folders.
% T columns: run | trial_idx | resp  (resp in {'resp1','resp2','timeout'})
% Adds [Run XX] prefix to the Data comment for GUI visibility.

% ---------- options ----------
p = inputParser;
addRequired(p,'T_or_path');
addRequired(p,'subjectName',@(s)ischar(s)||isstring(s));
addRequired(p,'runOffset',@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'SourceSuffix','_resp',@(s)ischar(s)||isstring(s));
addParameter(p,'TargetResp1','Resp1',@(s)ischar(s)||isstring(s));
addParameter(p,'TargetResp2','Resp2',@(s)ischar(s)||isstring(s));
addParameter(p,'DryRun',false,@islogical);
addParameter(p,'Verbose',true,@islogical);
addParameter(p,'MapMode','rank',@(s)any(strcmpi(s,{'rank','number'}))); % default: rank
addParameter(p,'PreviewN',10,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'TagRunInComment',true,@islogical);   % <-- NEW
parse(p,T_or_path,subjectName,runOffset,varargin{:});
opt = p.Results;
opt.MapMode = lower(string(opt.MapMode));

% ---------- load/normalize T ----------
if istable(opt.T_or_path)
    T = opt.T_or_path;
else
    T = readtable(char(opt.T_or_path));
end
assert(all(ismember({'run','trial_idx','resp'}, T.Properties.VariableNames)), ...
    'T must have columns: run, trial_idx, resp');
T.run  = string(T.run);
T.resp = lower(string(T.resp));
assert(all(ismember(unique(T.resp), ["resp1","resp2","timeout"])), ...
    'T.resp must be resp1/resp2/timeout');

% ---------- subject & target conditions ----------
[~, iSubject] = bst_get('Subject', char(opt.subjectName), 1); %#ok<ASGLU>
assert(~isempty(iSubject), 'Subject not found: %s', opt.subjectName);

iStudyR1 = ensure_condition(opt.subjectName, opt.TargetResp1);
iStudyR2 = ensure_condition(opt.subjectName, opt.TargetResp2);

destR1Rel = normalize_folder_rel(sprintf('data/%s/%s/', char(opt.subjectName), char(opt.TargetResp1)));
destR2Rel = normalize_folder_rel(sprintf('data/%s/%s/', char(opt.subjectName), char(opt.TargetResp2)));

% ---------- iterate runs ----------
runs_csv = unique(T.run, 'stable');
toCopy_R1 = {}; dst_R1 = {};
toCopy_R2 = {}; dst_R2 = {};
skipped = 0;

for k = 1:numel(runs_csv)
    runCsv = runs_csv(k);
    rnum   = local_parse_trailing_int(runCsv);              % "...run9" -> 9
    bsRun  = rnum + opt.runOffset;                          % apply offset
    condName   = sprintf('%s_Run_%02d%s', char(opt.subjectName), bsRun, char(opt.SourceSuffix));
    fullCondId = sprintf('%s/%s', char(opt.subjectName), condName);

    % Source condition
    [~, iStudy] = bst_get('StudyWithCondition', fullCondId);
    if isempty(iStudy)
        warning('Condition not found: %s', fullCondId);
        continue;
    end
    sStudy = bst_get('Study', iStudy);

    % Gather DB-relative file names
    files = {};
    if isfield(sStudy,'Data') && ~isempty(sStudy.Data)
        files = {sStudy.Data.FileName};
    end
    if isempty(files)
        warning('No data files in condition: %s', condName);
        continue;
    end
    files = normalize_dbrel_list(files);  % ensure 'data/' prefix and forward slashes

    % Keep only ...trialNNN.mat and sort by NNN
    [~, names, exts] = cellfun(@fileparts, files, 'UniformOutput', false);
    bases = strcat(names, exts);
    isData = ~cellfun('isempty', regexp(bases, 'trial(\d+)\.mat$', 'once'));
    files  = files(isData); bases = bases(isData);
    if isempty(files)
        warning('No "*trialNNN.mat" files in: %s', condName);
        continue;
    end
    trialNums = nan(numel(bases),1);
    for j = 1:numel(bases)
        tok = regexp(bases{j}, 'trial(\d+)\.mat$', 'tokens', 'once');
        if ~isempty(tok), trialNums(j) = str2double(tok{1}); end
    end
    keep = ~isnan(trialNums);
    files = files(keep); bases = bases(keep); trialNums = trialNums(keep);
    [trialNums, ord] = sort(trialNums); %#ok<ASGLU>
    files = files(ord); bases = bases(ord);

    % Select rows for this run (resp1/resp2 only)
    sel  = (T.run == runCsv) & ismember(T.resp, ["resp1","resp2"]);
    tIdx = T.trial_idx(sel);
    tRsp = T.resp(sel);
    if isempty(tIdx)
        if opt.Verbose
            fprintf('%s: No resp1/resp2 rows in T for this run (skipping).\n', condName);
        end
        continue;
    end

    % Pairing
    switch opt.MapMode
        case "rank"
            [~, ordT]  = sort(tIdx(:));                % rank within responded trials
            tRspSorted = tRsp(ordT);
            nFiles = numel(files); nRows = numel(tRspSorted);
            nPair  = min(nFiles, nRows);
            if nFiles ~= nRows
                warning('%s: %d files vs %d T rows; pairing by rank (min=%d).', ...
                        condName, nFiles, nRows, nPair);
            end
            for j = 1:nPair
                srcRel = files{j};                     % DB-relative
                [~, b, e] = fileparts(srcRel);
                newBase   = sanitize_name(sprintf('%s__%s%s', condName, b, e));  % includes Run_XX
                if tRspSorted(j) == "resp1"
                    dstRel = unique_relpath(destR1Rel, newBase, opt.DryRun);
                    toCopy_R1{end+1} = srcRel; %#ok<AGROW>
                    dst_R1{end+1}    = dstRel; %#ok<AGROW>
                else
                    dstRel = unique_relpath(destR2Rel, newBase, opt.DryRun);
                    toCopy_R2{end+1} = srcRel; %#ok<AGROW>
                    dst_R2{end+1}    = dstRel; %#ok<AGROW>
                end
            end
            skipped = skipped + abs(nFiles - nRows);

        case "number"
            % Map using original trial number embedded in filenames
            for j = 1:numel(tIdx)
                pos = find(trialNums == tIdx(j), 1, 'first');
                if isempty(pos)
                    skipped = skipped + 1;
                    continue;
                end
                srcRel = files{pos};
                [~, b, e] = fileparts(srcRel);
                newBase   = sanitize_name(sprintf('%s__%s%s', condName, b, e));
                if tRsp(j) == "resp1"
                    dstRel = unique_relpath(destR1Rel, newBase, opt.DryRun);
                    toCopy_R1{end+1} = srcRel; %#ok<AGROW>
                    dst_R1{end+1}    = dstRel; %#ok<AGROW>
                else
                    dstRel = unique_relpath(destR2Rel, newBase, opt.DryRun);
                    toCopy_R2{end+1} = srcRel; %#ok<AGROW>
                    dst_R2{end+1}    = dstRel; %#ok<AGROW>
                end
            end
    end
end

% ---------- plan / execute ----------
fprintf('[Plan] Copy %d -> %s, %d -> %s (skipped=%d)\n', ...
    numel(toCopy_R1), char(opt.TargetResp1), numel(toCopy_R2), char(opt.TargetResp2), skipped);

if opt.DryRun
    preview_moves('Resp1', toCopy_R1, dst_R1, opt.PreviewN);
    preview_moves('Resp2', toCopy_R2, dst_R2, opt.PreviewN);
    return;
end

copy_count = 0;
for i = 1:numel(toCopy_R1)
    srcAbs = rel2abs(toCopy_R1{i});
    dstAbs = rel2abs(dst_R1{i});
    [dstDir,~,~] = fileparts(dstAbs); if ~exist(dstDir,'dir'), mkdir(dstDir); end
    copyfile(srcAbs, dstAbs);
    copy_count = copy_count + 1;
end
for i = 1:numel(toCopy_R2)
    srcAbs = rel2abs(toCopy_R2{i});
    dstAbs = rel2abs(dst_R2{i});
    [dstDir,~,~] = fileparts(dstAbs); if ~exist(dstDir,'dir'), mkdir(dstDir); end
    copyfile(srcAbs, dstAbs);
    copy_count = copy_count + 1;
end


tag_run_prefix_abs(dst_R1);
tag_run_prefix_abs(dst_R2);
[s1,i1]=bst_get('StudyWithCondition','mcwa086_v1/Resp1');
[s2,i2]=bst_get('StudyWithCondition','mcwa086_v1/Resp2');
if ~isempty(i1), db_reload_studies(i1); end
if ~isempty(i2), db_reload_studies(i2); end
panel_protocols('UpdateTree');

% % ---------- add [Run XX] to comments (GUI) ----------
% if opt.TagRunInComment
%     tag_run_prefix(dst_R1);
%     tag_run_prefix(dst_R2);
% end

% Reload only the two aggregate studies and refresh GUI
if ~isempty(iStudyR1), db_reload_studies(iStudyR1); end
if ~isempty(iStudyR2), db_reload_studies(iStudyR2); end
panel_protocols('UpdateTree');

fprintf('Copied %d files. %d planned items were skipped.\n', copy_count, skipped);
end

% ---------- helpers ----------
function iStudy = ensure_condition(subjectName, cond)
fullCondId = sprintf('%s/%s', char(subjectName), char(cond));
[~, iStudy] = bst_get('StudyWithCondition', fullCondId);
if isempty(iStudy)
    [~, iSubject] = bst_get('Subject', char(subjectName), 1);
    iStudy = db_add_condition(iSubject, char(cond));
end
end

function n = local_parse_trailing_int(s)
m = regexp(char(s), '(\d+)$', 'tokens', 'once');
assert(~isempty(m), 'Cannot parse run number from "%s"', s);
n = str2double(m{1});
end

function relList = normalize_dbrel_list(relList)
relList = strrep(relList, '\','/');
for ii = 1:numel(relList)
    if ~(startsWith(relList{ii},'data/') || startsWith(relList{ii},'results/') || startsWith(relList{ii},'rawdata/'))
        relList{ii} = ['data/' relList{ii}];
    end
end
end

function folderRel = normalize_folder_rel(folderRel)
folderRel = strrep(folderRel,'\','/');
if ~(startsWith(folderRel,'data/') || startsWith(folderRel,'results/') || startsWith(folderRel,'rawdata/'))
    folderRel = ['data/' folderRel];
end
if ~endsWith(folderRel,'/'), folderRel = [folderRel '/']; end
end

function name = sanitize_name(name)
name = regexprep(name, '[^\w\.-]', '_');  % safe filename
end

function dstRel = unique_relpath(destFolderRel, baseName, isDryRun)
% Returns a DB-relative destination path that won't overwrite an existing file.
destFolderRel = normalize_folder_rel(destFolderRel);
candidateRel  = [destFolderRel baseName];

if isDryRun
    dstRel = candidateRel;   % no disk checks in dry-run
    return;
end

dstAbs = rel2abs(candidateRel);
[pth,b,ext] = fileparts(dstAbs);
k = 1;
finalAbs = dstAbs;
while exist(finalAbs,'file')
    finalAbs = fullfile(pth, sprintf('%s__c%02d%s', b, k, ext));
    k = k + 1;
end
dstRel = abs2rel(finalAbs);   % back to DB-relative
end

function absPath = rel2abs(relPath)
% Convert DB-relative ('data/...', 'results/...', 'rawdata/...') to absolute.
relPath = strrep(relPath,'\','/');
prot = bst_get('ProtocolInfo');
root = strrep(prot.STUDIES,'\','/'); % data root
if startsWith(relPath,'data/')
    tail = relPath(6:end);
elseif startsWith(relPath,'results/')
    tail = relPath(9:end);
elseif startsWith(relPath,'rawdata/')
    tail = relPath(9:end);
else
    tail = relPath;
end
absPath = fullfile(root, tail);
end

function relPath = abs2rel(absPath)
% Convert absolute path under Protocol STUDIES root back to 'data/...'
absPath = strrep(absPath,'\','/');
prot = bst_get('ProtocolInfo');
root = strrep(prot.STUDIES,'\','/');
assert(startsWith(absPath, root), 'Path outside STUDIES root.');
tail = absPath(numel(root)+2:end);  % strip "<root>/"
relPath = ['data/' tail];
relPath = strrep(relPath,'\','/');
end

function preview_moves(tag, srcList, dstList, N)
fprintf('  Examples -> %s (%d total):\n', tag, numel(srcList));
n = min(N, numel(srcList));
for i = 1:n
    fprintf('    %s  -->  %s\n', srcList{i}, dstList{i});
end
end

function tag_run_prefix(dstRelList)
% Add "[Run XX]" prefix to Data comments so run number is visible in GUI.
for i = 1:numel(dstRelList)
    fRel = dstRelList{i};                      % e.g., "data/mcwa086_v1/Resp1/xxx.mat" or "mcwa086_v1/Resp1/xxx.mat"
    tok  = regexp(fRel, 'Run_(\d+)', 'tokens', 'once');
    if isempty(tok), continue; end
    rstr = tok{1};                             % "02", "10", ...

    try
        % IMPORTANT: use absolute path for I/O to avoid "data/data" issues
        absF = rel2abs(fRel);

        D = in_bst_data(absF);                 % <-- absolute path
        oldC = ''; if isfield(D,'Comment') && ~isempty(D.Comment), oldC = D.Comment; end

        if contains(oldC, sprintf('[Run %s]', rstr))
            continue;                          % already tagged
        end
        baseC = regexprep(oldC, '^\[Run\s+\d+\]\s*', '');
        D.Comment = strtrim(sprintf('[Run %s] %s', rstr, strtrim(baseC)));

        bst_save(absF, D, 'v6');               % <-- absolute path
    catch ME
        warning('Could not tag comment for %s: %s', fRel, ME.message);
    end
end
end


function absPath = dbrel_to_abs(relOrAbs)
% Return an absolute path under the current protocol's STUDIES root.
relOrAbs = strrep(relOrAbs,'\','/');
prot = bst_get('ProtocolInfo');
root = strrep(prot.STUDIES,'\','/');  % e.g., /.../Pilot2_092025/data
if startsWith(relOrAbs, root)
    absPath = relOrAbs;                         % already absolute under STUDIES
elseif startsWith(relOrAbs, 'data/')
    absPath = fullfile(root, relOrAbs(6:end));  % strip "data/"
elseif startsWith(relOrAbs, 'results/')
    absPath = fullfile(root, relOrAbs(9:end));
elseif startsWith(relOrAbs, 'rawdata/')
    absPath = fullfile(root, relOrAbs(9:end));
else
    % treat as relative tail under STUDIES
    absPath = fullfile(root, relOrAbs);
end
absPath = strrep(absPath,'\','/');
end

function tag_run_prefix_abs(dstList)
% Update DataMat.Comment with "[Run XX]" using ABSOLUTE file paths only.
for i = 1:numel(dstList)
    fRel = strrep(dstList{i},'\','/');
    % Extract "XX" from "...Run_XX..."
    tok = regexp(fRel,'Run_(\d+)','tokens','once');
    if isempty(tok), continue; end
    rstr = tok{1};

    absF = dbrel_to_abs(fRel);
    if ~exist(absF,'file')
        warning('Missing on disk (skip): %s', absF);
        continue;
    end
    try
        S = load(absF,'-mat');                 % read Data file directly
        oldC = ''; if isfield(S,'Comment'), oldC = S.Comment; end
        if contains(oldC, sprintf('[Run %s]', rstr))
            continue;                           % already tagged
        end
        baseC = regexprep(oldC, '^\[Run\s+\d+\]\s*', '');
        newC  = strtrim(sprintf('[Run %s] %s', rstr, strtrim(baseC)));
        S.Comment = newC;

        % Save back; prefer bst_save, fall back to MATLAB save.
        try
            bst_save(absF, S, 'v6');
        catch
            save(absF, '-struct', 'S', '-v6');
        end
    catch ME
        warning('Tagging failed for %s: %s', absF, ME.message);
    end
end
end
