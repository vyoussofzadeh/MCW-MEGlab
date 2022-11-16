clear; clc, close('all'); warning off

% adding function path
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/BS_additions/BS_peakdetection/');

% adding BS path
bs_path = '/opt/matlab_toolboxes/Brainstorm3_2021/brainstorm3'; %BS 2021
addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

% peak roi detection
cfg = [];
cfg.npeaks = 3;
do_peakcoor_detection(cfg)