% Author: Vahab Youssof Zadeh
% Date (update): 03/30/2026
%
% Launcher for the MaxFilter + QC pipeline.
% Modes:
%   1) tSSS
%   2) tSSS + trans
%   3) tSSS + trans + movecomp

clc; clear; close all;

addpath('/MEG_data/MCW_pipeline/Preprocess/maxfilter_for_reseach_data/func/')

% -------------------------------------------------------------------------
% Enter session folder
% -------------------------------------------------------------------------
disp('---');
disp('Enter the session/data folder:');
disp('Example: /MEG_data/neuro-data/tacs_stroke/test_test/260202');
sessionDir = strtrim(input('Session folder: ', 's'));

if isempty(sessionDir)
    error('Session folder was not provided.');
end

if exist(sessionDir, 'dir') ~= 7
    error('Session folder not found: %s', sessionDir);
end

origDir = fullfile(sessionDir, 'orig');
if exist(origDir, 'dir') ~= 7
    error('orig folder not found: %s', origDir);
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
        bashScript = '/MEG_data/MCW_pipeline/Preprocess/maxfilter_for_reseach_data/vectorview/vview_preProc_maxfilter.sh';
        tsssDir    = fullfile(sessionDir, 'tsss');

    case 2
        bashScript = '/MEG_data/MCW_pipeline/Preprocess/maxfilter_for_reseach_data/vectorview/vview_preProc_maxfilter_trans.sh';
        tsssDir    = fullfile(sessionDir, 'tsss_tans');

        disp('---');
        disp('Enter the baseline FIF file used for transformation (-trans):');
        disp('Example: mcwa117_pSTM_run1_raw.fif');
        baselineFile = strtrim(input('Baseline FIF file: ', 's'));

        if isempty(baselineFile)
            error('Baseline FIF file was not provided.');
        end
    otherwise
        error('Invalid mode. Please choose 1, 2, or 3.');
end

% -------------------------------------------------------------------------
% Checks
% -------------------------------------------------------------------------
if exist(bashScript, 'file') ~= 2
    error('Bash script not found: %s', bashScript);
end

if ~isempty(baselineFile)
    baselinePath = fullfile(origDir, baselineFile);
    if exist(baselinePath, 'file') ~= 2
        error('Baseline FIF file not found in orig/: %s', baselinePath);
    end
else
    baselinePath = '';
end

logDir = fullfile(sessionDir, 'pipeline_logs');

fprintf('\nStarting MaxFilter pipeline...\n');
fprintf('Session folder : %s\n', sessionDir);
fprintf('orig folder    : %s\n', origDir);
fprintf('Mode           : %d\n', modeID);
fprintf('Output folder  : %s\n', tsssDir);
fprintf('Bash script    : %s\n', bashScript);

if ~isempty(baselinePath)
    fprintf('Baseline file  : %s\n', baselinePath);
else
    fprintf('Baseline file  : not required\n');
end
fprintf('\n');

% -------------------------------------------------------------------------
% Run pipeline
% -------------------------------------------------------------------------
Summary = runmaxfilter(sessionDir, baselinePath, bashScript, tsssDir, logDir);

fprintf('\nPipeline finished.\n');
disp('Summary output:');
disp(Summary);

