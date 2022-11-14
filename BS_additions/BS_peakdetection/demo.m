clear; clc, close('all'); warning off


restoredefaultpath
% clear
clc
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/BS_additions/BS_peakdetection/');

bs_path = '/opt/matlab_toolboxes/brainstorm3';
addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

cfg = [];
cfg.npeaks = 8;
cfg.BSdir = 'Group_analysis/Five_intervals';
do_peakcoor_detection(cfg)