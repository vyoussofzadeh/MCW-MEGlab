clear, clc
%% === CONFIG ===
rootDir    = '/MEG_data/AHW_SpikeAnalysis';
sheetDir   = fullfile(rootDir, 'Review_sheet');
dataDir    = fullfile(rootDir, 'MEG_data/sss');     % where .fif live
destRoot   = fullfile(rootDir, 'Reviewer');         % Adi/Josh/Manoj/Pradeep live here

useSymlinks       = false;   % false = copy actual files, true = ln -s (saves space)
includeEventFiles = false;    % true  = include "*-eve.fif" alongside raw .fif
byRunSubfolder    = false;   % true  = Reviewer/<Reviewer>/<Subject>/RunXX/
overwriteExisting = true;    % true  = replace existing file at destination

%% === Load the newest assignments CSV ===
d = dir(fullfile(sheetDir, 'review_assignments*.csv'));
assert(~isempty(d), 'No assignment CSV found in %s', sheetDir);
[~,iNewest] = max([d.datenum]);
assignCsv   = fullfile(d(iNewest).folder, d(iNewest).name);
fprintf('Using assignments: %s\n', assignCsv);

T = readtable(assignCsv, 'TextType','string');

% Expected columns: Subject, (RunNum or RunLabel), FileName/FullPath (optional), Physician1, Physician2
% Ensure we have RunNum
if ~ismember('RunNum', T.Properties.VariableNames)
    if ismember('RunLabel', T.Properties.VariableNames)
        T.RunNum = double(extractAfter(T.RunLabel, "Run"));
    else
        % Parse from FileName or FullPath as a fallback
        srcCol = [];
        if ismember('FileName', T.Properties.VariableNames), srcCol = T.FileName; end
        if isempty(srcCol) && ismember('FullPath', T.Properties.VariableNames), srcCol = string(T.FullPath); end
        assert(~isempty(srcCol), 'Need RunNum or a parsable FileName/FullPath in the assignment CSV.');
        parsed = regexp(srcCol, '_(?:run|Run)0?(\d+)_', 'tokens', 'once');
        T.RunNum = nan(height(T),1);
        for k = 1:height(T)
            if ~isempty(parsed{k}), T.RunNum(k) = str2double(parsed{k}{1}); end
        end
        assert(all(~isnan(T.RunNum)), 'Could not parse RunNum for all rows.');
    end
end

% Long format: one row per reviewer per run
R1 = T(:, {'Subject','RunNum'}); R1.Reviewer = T.Physician1;
R2 = T(:, {'Subject','RunNum'}); R2.Reviewer = T.Physician2;
Rows = [R1; R2];

% Drop any incomplete rows (just in case)
Rows = Rows(~(Rows.Subject=="" | isnan(Rows.RunNum) | Rows.Reviewer==""), :);

% Unique jobs so we don't duplicate copies per reviewer/run
[~, uix] = unique(Rows(:, {'Reviewer','Subject','RunNum'}), 'rows', 'stable');
Jobs = Rows(uix, :);

% Ensure reviewer folders exist
reviewers = unique(Jobs.Reviewer);
for r = 1:numel(reviewers)
    rdir = fullfile(destRoot, reviewers(r));
    if ~exist(rdir, 'dir'), mkdir(rdir); end
end

%% === Helper to find matching .fif files for a subject+run ===
% Matches both Run06 and Run6, case-insensitive on "Run"
findRunFiles = @(subj, runNum) local_find_run_files(dataDir, subj, runNum, includeEventFiles);

%% === Copy/link files ===
nOK = 0; nMiss = 0; nFail = 0;
logRows = strings(0,6);  % Reviewer,Subject,Run,Src,Dest,Action

for k = 1:height(Jobs)
    subj  = string(Jobs.Subject(k));
    runN  = Jobs.RunNum(k);
    rev   = string(Jobs.Reviewer(k));

    % where to drop them
%     dstBase = fullfile(destRoot, rev, subj);
    dstBase = fullfile(destRoot, rev, 'sss/MEG_data');
    if byRunSubfolder
        dstBase = fullfile(dstBase, sprintf('Run%02d', runN));
    end
    if ~exist(dstBase, 'dir'), mkdir(dstBase); end

    % locate .fif files for this subject/run
    files = findRunFiles(subj, runN);
    if isempty(files)
        warning('No .fif found for %s Run%d', subj, runN);
        nMiss = nMiss + 1;
        logRows(end+1,:) = [rev, subj, string(runN), "", "", "MISSING"]; %#ok<SAGROW>
        continue;
    end

    % place all matches
    for f = 1:numel(files)
        src = fullfile(dataDir, files(f).name);
        dst = fullfile(dstBase, files(f).name);

        if useSymlinks
            % create/update symlink
            if overwriteExisting && exist(dst, 'file'), delete(dst); end
            cmd = sprintf('ln -s "%s" "%s"', src, dst);
            if exist(dst, 'file')
                % If a stale file exists and overwrite=false, skip
                action = "SKIP";
                if ~overwriteExisting
                    logRows(end+1,:) = [rev, subj, string(runN), src, dst, action]; % #ok<SAGROW>
                    continue;
                end
            end
            [status, msg] = system(cmd);
            if status ~= 0
                warning('Symlink failed: %s -> %s (%s). Falling back to copy.', src, dst, strtrim(msg));
                [ok, msg2] = copyfile(src, dst, 'f');
                if ok
                    nOK = nOK + 1; action = "COPY(FB)";
                else
                    warning('Copy failed: %s', msg2); nFail = nFail + 1; action = "FAIL";
                end
            else
                nOK = nOK + 1; action = "LINK";
            end
        else
            % real copy
            [ok,msg2] = copyfile_force(src, dst, overwriteExisting);

%             [ok, msg2] = copyfile(src, dst, overwriteExisting * 'f');
            if ok
                nOK = nOK + 1; action = "COPY";
            else
                warning('Copy failed: %s -> %s (%s)', src, dst, msg2);
                nFail = nFail + 1; action = "FAIL";
            end
        end

        logRows(end+1,:) = [rev, subj, string(runN), src, dst, action]; %#ok<SAGROW>
    end
end

%% === Summary + optional log CSV ===
fprintf('Placed %d items (%d missing subject/run, %d failed).\n', nOK, nMiss, nFail);

logT = array2table(logRows, 'VariableNames', ...
    {'Reviewer','Subject','Run','Source','Destination','Action'});
logCsv = fullfile(sheetDir, 'copy_log.csv');
writetable(logT, logCsv);
fprintf('Log written to %s\n', logCsv);

%% === Local function(s) ===
function files = local_find_run_files(dataDir, subject, runNum, includeEvent)
% Find .fif for subject+run (supports Run6 and Run06), optionally include *-eve.fif
% Returns an array of dir() results (may be empty)

    % Build two patterns (zero-padded and non-padded)
    pat1 = sprintf('%s_Run%02d_*.fif', subject, runNum);
    pat2 = sprintf('%s_Run%d_*.fif',  subject, runNum);

    L1 = dir(fullfile(dataDir, pat1));
    L2 = dir(fullfile(dataDir, pat2));

    files = [L1; L2];

    % Deduplicate by name
    if ~isempty(files)
        [~, ia] = unique({files.name}, 'stable');
        files = files(ia);
    end

    if ~includeEvent && ~isempty(files)
        keep = ~endsWith({files.name}, '-eve.fif');
        files = files(keep);
    end
end

function [ok,msg] = copyfile_force(src,dst,force)
    if force, [ok,msg] = copyfile(src,dst,'f'); else, [ok,msg] = copyfile(src,dst); end
end

