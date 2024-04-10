
%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
addpath('./run')
addpath('./function')
Run_setpath

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other/')

%%
LI_analysis_label = {'DICS_indirect','DICS_directcontrast','LCMV_anim_vs_Symb','-','DICS_anim', 'DICS_contrast_prestim', 'dSPM_contrast'};


for i = 1:length(LI_analysis_label)
    disp([num2str(i) ') ' LI_analysis_label{i}]);
end

LI_analysis = input('');

%%
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/';
cd(data_save_dir)

%% Run Brainstorm
Run_BS
brainstorm

%%
Run_load_surface_template

%% Read fMRI lats
fmri_LIs = ecpfunc_read_fmri_lat();

%%
datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';
BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data';

cfg = [];
cfg.protocol = protocol;
cfg.datadir = datadir;
cfg.BS_data_dir = BS_data_dir;

switch LI_analysis
    case {1,5}
        cfg.datatag = 'wDICS_baseline_18_4';
        S_data = ecpfunc_read_sourcemaps_dics(cfg);
    case 2
        cfg.datatag = 'wDICS_contrast_18_4';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 3
        cfg.datamask = fullfile('./Group_analysis/LCMV/results_average*.mat');
        S_data = ecpfunc_read_sourcemaps(cfg);
    case 6
        cfg.datatag = 'PSTwDICS_contrast_18_4';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 7
        cfg.datatag = 'dSPM_contrast';
        S_data = ecpfunc_read_sourcemaps_contrast(cfg);
end

%% TLE side (PT only)
patn_neuropsych_data = ecpfunc_read_patn_neuropsych();

%% Subject demog details
switch LI_analysis
    case {1,3 ,5}
        cfg = []; cfg.subjs_3 = S_data.subjs_3; cfg.subjs_2 = S_data.subjs_2;
        cfg.sFiles_3 = S_data.sFiles_3; cfg.sFiles_2 = S_data.sFiles_2;
        sub_demog_data = ecpfunc_read_sub_demog(cfg);
    case {2, 6, 7}
        cfg = []; cfg.subjs = S_data.subjs;
        cfg.sFiles = S_data.sFiles_32;
        sub_demog_data = ecpfunc_read_sub_demog_contrast(cfg);
end

%% Inter-subject (group) averaging,
switch LI_analysis
    case {1,3, 5}
        disp('1: Anim, Ctrl')
        disp('2: Anim, Patn')
        disp('3: Symbol, Ctrl')
        disp('4: Symbol, Patn')
        cfg = [];
        cfg.sub_demog_data = sub_demog_data;
        cfg.select_data = 1; S_data_anim_hc = ecpfunc_select_data(cfg);
        cfg.select_data = 2; S_data_anim_pt = ecpfunc_select_data(cfg);
        cfg.select_data = 3; S_data_symb_hc = ecpfunc_select_data(cfg);
        cfg.select_data = 4; S_data_symb_pt = ecpfunc_select_data(cfg);
        [LI_pt_ID,~,~] = intersect(S_data_anim_pt.sFiles_subid, S_data_symb_pt.sFiles_subid);
        [LI_hc_ID,~,~] = intersect(S_data_anim_hc.sFiles_subid, S_data_symb_hc.sFiles_subid);        
        
    case {2,4,6}
        disp('1: Ctrl')
        disp('2: Patn')
        cfg = [];
        cfg.sub_demog_data = sub_demog_data;
        cfg.select_data = 1; S_data_hc = ecpfunc_select_contrast_data(cfg);
        cfg.select_data = 2; S_data_pt = ecpfunc_select_contrast_data(cfg);
        LI_pt_ID = S_data_pt.sFiles_subid;
end

%%
switch LI_analysis
    case 1
        
    case 5
        sFiles_in = S_data_anim_pt.sFiles_in;
    case {2,6}
        sFiles_in = S_data_pt.sFiles_in;
end

%% avg. PT
bst_process('CallProcess', 'process_average', sFiles_in, [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'scalenormalized', 0);

% brainstormData = load(sFiles.FileName);

%% Stats 
% Process: t-test zero [-0.300s,2.100s]          H0:(X=0), H1:(X<>0)
bst_process('CallProcess', 'process_test_parametric1', sFiles_in, [], ...
    'timewindow',    [-0.3, 2.1], ...
    'scoutsel',      {}, ...
    'scoutfunc',     1, ...  % Mean
    'isnorm',        0, ...
    'avgtime',       0, ...
    'Comment',       '', ...
    'test_type',     'ttest_onesample', ...  % One-sample Student's t-test    X~N(m,s)t = mean(X) ./ std(X) .* sqrt(n)      df=n-1
    'tail',          'two');  % Two-tailed

% brainstormData = load(sFiles.FileName);

%%
sFilesB = db_template('importfile');

% Process: FT t-test unequal cluster [0.000s,2.100s]          H0:(A=B), H1:(A<>B)
bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_in, sFilesB, ...
    'timewindow',     [7.21644966e-16, 2.1], ...
    'scoutsel',       {}, ...
    'scoutfunc',      1, ...  % Mean
    'isabs',          0, ...
    'avgtime',        0, ...
    'randomizations', 1000, ...
    'statistictype',  1, ...  % Independent t-test
    'tail',           'two', ...  % Two-tailed
    'correctiontype', 2, ...  % cluster
    'minnbchan',      0, ...
    'clusteralpha',   0.05);

%% Interval Avg.
disp('import avg/stats data brainstormData Matlab variable')
pause,

% Assuming your data is stored in a variable named 'brainstormData'
startTime = -0.2; % Example start time in seconds
endTime = 0;   % Example end time in seconds

startTime = 0; % Example start time in seconds
endTime = .200;   % Example end time in seconds

startTime = .200; % Example start time in seconds
endTime = .400;   % Example end time in seconds

startTime = .400; % Example start time in seconds
endTime = .600;   % Example end time in seconds

startTime = .600; % Example start time in seconds
endTime = .800;   % Example end time in seconds

startTime = .800; % Example start time in seconds
endTime = 1;   % Example end time in seconds


startTime = -0.2; % Example start time in seconds
endTime = 0.1;   % Example end time in seconds

startTime = 0.1; % Example start time in seconds
endTime = 0.4;   % Example end time in seconds

startTime = .400; % Example start time in seconds
endTime = .700;   % Example end time in seconds

startTime = .700; % Example start time in seconds
endTime = 1;   % Example end time in seconds

startTime = 1; % Example start time in seconds
endTime = 1.3;   % Example end time in seconds

startTime = 1.3; % Example start time in seconds
endTime = 1.6;   % Example end time in seconds

startTime = 1.6; % Example start time in seconds
endTime = 1.9;   % Example end time in seconds

newbrainstormData = calculateIntervalAverage(brainstormData, startTime, endTime);
disp('import newbrainstormData using BS GUI')
