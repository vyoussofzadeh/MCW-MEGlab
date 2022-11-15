clear; clc, close('all'); warning off

clc
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/BS_additions/BS_peakdetection/');

bs_path = '/opt/matlab_toolboxes/Brainstorm3_2021/brainstorm3'; %BS 2021

addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

cfg = [];
cfg.npeaks = 3;
do_peakcoor_detection(cfg)