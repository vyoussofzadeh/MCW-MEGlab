% Author: Vahab Youssof Zadeh
% Date: 03/23/2026
%
% Simple launcher for the MaxFilter + QC pipeline.
% The script asks for:
%   1) session/data directory
%   2) baseline FIF file used for -trans
% Then it calls:
%   Summary = runmaxfilter(sessionDir, baselineFile, bashScript);

clc; clear; close all;

% Path to the bash MaxFilter wrapper
bashScript = '/MEG_data/MCW_pipeline/Preprocess/maxfilter_for_reseach_data/preProc_maxfilter.sh';

% Example:
% sessionDir   = /MEG_data/neuro-data/tacs_stroke/test_test/260202
% baselineFile = mcwa117_pSTM_run1_raw.fif

disp('Enter the session/data folder:');
disp('Example: /MEG_data/neuro-data/tacs_stroke/test_test/260202');
sessionDir = strtrim(input('Session folder: ', 's'));

if isempty(sessionDir)
    error('Session folder was not provided.');
end

if exist(sessionDir, 'dir') ~= 7
    error('Session folder not found: %s', sessionDir);
end

disp('Enter the baseline FIF file used for transformation (-trans):');
disp('Example: mcwa117_pSTM_run1_raw.fif');
baselineFile = strtrim(input('Baseline FIF file: ', 's'));

if isempty(baselineFile)
    error('Baseline FIF file was not provided.');
end

if exist(bashScript, 'file') ~= 2
    error('Bash script not found: %s', bashScript);
end

fprintf('\nStarting MaxFilter pipeline...\n');
fprintf('Session folder : %s\n', sessionDir);
fprintf('Baseline file  : %s\n', baselineFile);
fprintf('Bash script    : %s\n\n', bashScript);

tsssDir = fullfile(sessionDir, 'tsss_tans_mvcomp');
logDir  = fullfile(sessionDir, 'pipeline_logs');

Summary = runmaxfilter(sessionDir, baselineFile, bashScript, tsssDir, logDir);

fprintf('\nPipeline finished.\n');
disp('Summary output:');
disp(Summary);