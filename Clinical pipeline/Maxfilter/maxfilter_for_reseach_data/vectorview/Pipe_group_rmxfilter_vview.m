% Author: Vahab Youssof Zadeh
% Date (update): 03/31/2026
%
% Batch launcher for the MaxFilter + QC pipeline.
% Modes:
%   1) tSSS
%   2) tSSS + trans
%
% Assumptions:
% - Each session folder contains an orig/ folder with raw FIF files
% - runmaxfilter(sessionDir, baselinePath, bashScript, tsssDir, logDir) exists
% - For mode 2, the baseline filename is the same across sessions

clc; clear; close all;

addpath('/MEG_data/MCW_pipeline/Preprocess/maxfilter_for_reseach_data/func/')

% -------------------------------------------------------------------------
% Enter study root folder
% -------------------------------------------------------------------------
disp('---');
disp('Enter the ROOT folder that contains all subject/session folders:');
disp('Example: /MEG_data/neuro-data/tacs_stroke');
rootDir = strtrim(input('Root folder: ', 's'));

if isempty(rootDir)
    error('Root folder was not provided.');
end

if exist(rootDir, 'dir') ~= 7
    error('Root folder not found: %s', rootDir);
end

% -------------------------------------------------------------------------
% Select mode
% -------------------------------------------------------------------------
disp('---');
disp('Select preprocessing mode:');
disp('  1 = tSSS');
disp('  2 = tSSS + trans');
modeID = input('Mode: ');

baselineFile = '';

switch modeID
    case 1
        bashScript    = '/MEG_data/MCW_pipeline/Preprocess/maxfilter_for_reseach_data/vectorview/vview_preProc_maxfilter.sh';
        outFolderName = 'tsss';
        
    case 2
        bashScript    = '/MEG_data/MCW_pipeline/Preprocess/maxfilter_for_reseach_data/vectorview/vview_preProc_maxfilter_trans.sh';
        outFolderName = 'tsss_tans';
        
    otherwise
        error('Invalid mode. Please choose 1 or 2.');
end

if exist(bashScript, 'file') ~= 2
    error('Bash script not found: %s', bashScript);
end

% -------------------------------------------------------------------------
% Find all session folders
% -------------------------------------------------------------------------
sessionDirs = find_session_dirs(rootDir);

% Keep only HC sessions, exclude MCW/patient sessions
keepMask = false(size(sessionDirs));

for i = 1:numel(sessionDirs)
    thisPath = lower(sessionDirs{i});

    isHC  = ~isempty(regexp(thisPath, '[/\\]hc[^/\\]*[/\\]', 'once'));
    isMCW = ~isempty(regexp(thisPath, '[/\\]mcw[^/\\]*[/\\]', 'once'));

    keepMask(i) = isHC && ~isMCW;
end

sessionDirs = sessionDirs(keepMask);

if isempty(sessionDirs)
    error('No session folders with orig/ were found under: %s', rootDir);
end

fprintf('\nFound %d session(s) with orig/ folders.\n', numel(sessionDirs));

% -------------------------------------------------------------------------
% Process all sessions
% -------------------------------------------------------------------------
BatchSummary = struct( ...
    'sessionDir', {}, ...
    'status', {}, ...
    'modeID', {}, ...
    'nRaw', {}, ...
    'nDoneBefore', {}, ...
    'nMissingBefore', {}, ...
    'baselinePath', {}, ...
    'message', {});

if modeID == 2
    disp('---');
    disp('Choose baseline selection method:');
    disp('  1 = first collected (task) raw run in each session');
    disp('  2 = specific filename in each session');
    disp('  3 = choose manually per session');
    baselineMode = input('Baseline mode: ');
    
    baselinePattern = '';
    if baselineMode == 2
        baselinePattern = strtrim(input('Baseline filename: ', 's'));
        if isempty(baselinePattern)
            error('Baseline filename was not provided.');
        end
    elseif ~ismember(baselineMode, [1 2 3])
        error('Invalid baseline mode. Choose 1, 2, or 3.');
    end
else
    baselineMode = 0;
    baselinePattern = '';
end

for iSess = 18:numel(sessionDirs)
    
    disp(sessionDirs(iSess))
    
    sessionDir = sessionDirs{iSess};
    origDir    = fullfile(sessionDir, 'orig');
    tsssDir    = fullfile(sessionDir, outFolderName);
    logDir     = fullfile(sessionDir, 'pipeline_logs');
    
    cd(sessionDir)
    
    fprintf('\n=====================================================\n');
    fprintf('Session %d / %d\n', iSess, numel(sessionDirs));
    fprintf('Session folder: %s\n', sessionDir);
    
    rawFiles = dir(fullfile(origDir, '*_raw.fif'));
    rawFiles = rawFiles(~[rawFiles.isdir]);
    
    if isempty(rawFiles)
        fprintf('No raw FIF files found. Skipping.\n');
        BatchSummary(end+1) = make_batch_row(sessionDir, 'skip_no_raw', modeID, 0, 0, 0, '', 'No raw FIF files found'); %#ok<SAGROW>
        continue;
    end
    
    [nDoneBefore, nMissingBefore] = count_done_vs_missing(origDir, tsssDir);
    
    fprintf('Raw runs found      : %d\n', numel(rawFiles));
    fprintf('Already processed   : %d\n', nDoneBefore);
    fprintf('Still missing       : %d\n', nMissingBefore);
    
    % Skip if everything is already done
    if nMissingBefore == 0
        fprintf('All outputs already exist. Skipping session.\n');
        BatchSummary(end+1) = make_batch_row(sessionDir, 'skip_done', modeID, numel(rawFiles), nDoneBefore, nMissingBefore, '', 'All outputs already exist'); %#ok<SAGROW>
        continue;
    end
    
    
    if modeID == 2
        [baselinePath, msg] = choose_session_baseline(origDir, rawFiles, baselineMode, baselinePattern);
        fprintf('%s\n', msg);
        
        if isempty(baselinePath)
            BatchSummary(end+1) = make_batch_row( ...
                sessionDir, 'skip_no_baseline', modeID, numel(rawFiles), ...
                nDoneBefore, n1MissingBefore, '', msg); %#ok<SAGROW>
            continue;
        end
    else
        baselinePath = '';
    end
    
    fprintf('Output folder       : %s\n', tsssDir);
    fprintf('Bash script         : %s\n', bashScript);
    if ~isempty(baselinePath)
        fprintf('Baseline file       : %s\n', baselinePath);
    else
        fprintf('Baseline file       : not required\n');
    end
    
    try
        Summary = runmaxfilter(sessionDir, baselinePath, bashScript, tsssDir, logDir);
        
        [nDoneAfter, nMissingAfter] = count_done_vs_missing(origDir, tsssDir);
        
        if nMissingAfter == 0
            statusTxt = 'done';
            msgTxt = 'Session completed successfully';
        elseif nDoneAfter > nDoneBefore
            statusTxt = 'partial_progress';
            msgTxt = sprintf('Some runs processed, but %d still missing', nMissingAfter);
        else
            statusTxt = 'no_progress';
            msgTxt = 'runmaxfilter returned but no additional outputs were detected';
        end
        
        BatchSummary(end+1) = make_batch_row(sessionDir, statusTxt, modeID, numel(rawFiles), nDoneBefore, nMissingBefore, baselinePath, msgTxt); %#ok<SAGROW>
        
    catch ME
        fprintf('ERROR: %s\n', ME.message);
        BatchSummary(end+1) = make_batch_row(sessionDir, 'error', modeID, numel(rawFiles), nDoneBefore, nMissingBefore, baselinePath, ME.message); %#ok<SAGROW>
    end
end

% -------------------------------------------------------------------------
% Save batch summary
% -------------------------------------------------------------------------
batchMat = fullfile(rootDir, sprintf('batch_maxfilter_mode%d.mat', modeID));
save(batchMat, 'BatchSummary');

batchCsv = fullfile(rootDir, sprintf('batch_maxfilter_mode%d.csv', modeID));
write_batch_summary_csv(batchCsv, BatchSummary);

fprintf('\n=====================================================\n');
fprintf('Batch processing finished.\n');
fprintf('Summary MAT: %s\n', batchMat);
fprintf('Summary CSV: %s\n', batchCsv);

disp(struct2table(BatchSummary));

function sessionDirs = find_session_dirs(rootDir)

sessionDirs = {};

if exist(fullfile(rootDir, 'orig'), 'dir') == 7
    sessionDirs{end+1} = rootDir;
end

d = dir(rootDir);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.','..'}));

for i = 1:numel(d)
    childDir = fullfile(rootDir, d(i).name);
    childSessions = find_session_dirs(childDir);
    sessionDirs = [sessionDirs, childSessions]; %#ok<AGROW>
end
end

function [nDone, nMissing] = count_done_vs_missing(origDir, outDir)

rawFiles = dir(fullfile(origDir, '*_raw.fif'));
rawFiles = rawFiles(~[rawFiles.isdir]);

nDone = 0;
nMissing = 0;

for i = 1:numel(rawFiles)
    outFile = fullfile(outDir, rawFiles(i).name);

    isDone = false;
    if exist(outFile, 'file') == 2
        d = dir(outFile);
        if ~isempty(d) && d(1).bytes > 0
            isDone = true;
        end
    end

    if isDone
        nDone = nDone + 1;
    else
        nMissing = nMissing + 1;
    end
end
end

function S = make_batch_row(sessionDir, statusTxt, modeID, nRaw, nDoneBefore, nMissingBefore, baselinePath, msgTxt)

S = struct( ...
    'sessionDir', sessionDir, ...
    'status', statusTxt, ...
    'modeID', modeID, ...
    'nRaw', nRaw, ...
    'nDoneBefore', nDoneBefore, ...
    'nMissingBefore', nMissingBefore, ...
    'baselinePath', baselinePath, ...
    'message', msgTxt);
end

function write_batch_summary_csv(csvFile, BatchSummary)

fid = fopen(csvFile, 'wt');
if fid == -1
    error('Cannot write CSV file: %s', csvFile);
end
cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid, 'sessionDir,status,modeID,nRaw,nDoneBefore,nMissingBefore,baselinePath,message\n');

for i = 1:numel(BatchSummary)
    s = BatchSummary(i);
    fprintf(fid, '%s,%s,%d,%d,%d,%d,%s,%s\n', ...
        csv_escape_local(s.sessionDir), ...
        csv_escape_local(s.status), ...
        s.modeID, ...
        s.nRaw, ...
        s.nDoneBefore, ...
        s.nMissingBefore, ...
        csv_escape_local(s.baselinePath), ...
        csv_escape_local(s.message));
end
end

function out = csv_escape_local(str)
if isempty(str)
    out = '""';
    return;
end
str = strrep(str, '"', '""');
out = ['"' str '"'];
end

function [baselinePath, msg] = choose_session_baseline(origDir, rawFiles, baselineMode, baselinePattern)

baselinePath = '';
msg = '';

% Exclude these from baseline candidates (case-insensitive)
excludeKeywords = {'emptyroom', 'Noise'};

% Build keep mask
keepMask = true(size(rawFiles));
for i = 1:numel(rawFiles)
    thisName = rawFiles(i).name;
    for k = 1:numel(excludeKeywords)
        if ~isempty(regexpi(thisName, excludeKeywords{k}))   % case-insensitive
            keepMask(i) = false;
            break;
        end
    end
end

rawFilesSel = rawFiles(keepMask);

if isempty(rawFilesSel)
    msg = 'No valid baseline candidates left after excluding emptyroom and ERNoise.';
    return;
end

switch baselineMode
    case 1
        [~, idx] = sort([rawFilesSel.datenum], 'ascend');
        rawFilesSel = rawFilesSel(idx);

        baselinePath = fullfile(origDir, rawFilesSel(1).name);
        msg = sprintf('Using first collected run as baseline: %s', baselinePath);

    case 2
        isExcluded = false;
        for k = 1:numel(excludeKeywords)
            if ~isempty(regexpi(baselinePattern, excludeKeywords{k}))   % case-insensitive
                isExcluded = true;
                break;
            end
        end

        if isExcluded
            msg = sprintf('Requested baseline is excluded by rule: %s', baselinePattern);
            return;
        end

        candidate = fullfile(origDir, baselinePattern);
        if exist(candidate, 'file') == 2
            baselinePath = candidate;
            msg = sprintf('Using requested baseline: %s', baselinePath);
        else
            msg = sprintf('Requested baseline not found: %s', candidate);
        end

    case 3
        [~, idx] = sort([rawFilesSel.datenum], 'ascend');
        rawFilesSel = rawFilesSel(idx);

        fprintf('\nAvailable raw runs in this session:\n');
        for i = 1:numel(rawFilesSel)
            fprintf('  %d = %s   [%s]\n', i, rawFilesSel(i).name, datestr(rawFilesSel(i).datenum));
        end

        idxChoice = input('Choose baseline run number (empty to skip session): ');

        if isempty(idxChoice)
            msg = 'No baseline selected.';
            return;
        end

        if ~isscalar(idxChoice) || idxChoice < 1 || idxChoice > numel(rawFilesSel) || idxChoice ~= round(idxChoice)
            msg = 'Invalid baseline selection.';
            return;
        end

        baselinePath = fullfile(origDir, rawFilesSel(idxChoice).name);
        msg = sprintf('Using selected baseline: %s', baselinePath);

    otherwise
        msg = 'Invalid baseline mode.';
end

% function [baselinePath, msg] = choose_session_baseline(origDir, rawFiles, baselineMode, baselinePattern)
% 
% baselinePath = '';
% msg = '';
% 
% switch baselineMode
%     case 1
%         % First collected raw run based on file timestamp
%         [~, idx] = sort([rawFiles.datenum], 'ascend');
%         rawFiles = rawFiles(idx);
%         
%         baselinePath = fullfile(origDir, rawFiles(1).name);
%         msg = sprintf('Using first collected run as baseline: %s', baselinePath);
%         
%     case 2
%         candidate = fullfile(origDir, baselinePattern);
%         if exist(candidate, 'file') == 2
%             baselinePath = candidate;
%             msg = sprintf('Using requested baseline: %s', baselinePath);
%         else
%             msg = sprintf('Requested baseline not found: %s', candidate);
%         end
%         
%     case 3
%         % Manual per session, shown in collection order
%         [~, idx] = sort([rawFiles.datenum], 'ascend');
%         rawFiles = rawFiles(idx);
%         
%         fprintf('\nAvailable raw runs in this session:\n');
%         for i = 1:numel(rawFiles)
%             fprintf('  %d = %s   [%s]\n', i, rawFiles(i).name, datestr(rawFiles(i).datenum));
%         end
%         
%         idxChoice = input('Choose baseline run number (empty to skip session): ');
%         
%         if isempty(idxChoice)
%             msg = 'No baseline selected.';
%             return;
%         end
%         
%         if ~isscalar(idxChoice) || idxChoice < 1 || idxChoice > numel(rawFiles) || idxChoice ~= round(idxChoice)
%             msg = 'Invalid baseline selection.';
%             return;
%         end
%         
%         baselinePath = fullfile(origDir, rawFiles(idxChoice).name);
%         msg = sprintf('Using selected baseline: %s', baselinePath);
%         
%     otherwise
%         msg = 'Invalid baseline mode.';
% end
end