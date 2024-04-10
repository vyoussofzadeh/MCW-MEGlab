
%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 08/09/2023

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
LI_method_label = {'Magnitude', 'Counting','Bootstrapping'};

disp('1: Magnitude')
disp('2: Counting')
disp('3: Bootstrapping')
LI_method = input(':');

%%
data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results';
cd(data_save_dir)

%%
save_dir = fullfile(data_save_dir,LI_analysis_label{LI_analysis}, LI_method_label{LI_method});
checkOrCreateDir(save_dir)
cd(save_dir)

LI_hc = load('LI_Ctrl');
LI_pt = load('LI_Patn');

%% Run Brainstorm
Run_BS

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

%% Subject demog details
switch LI_analysis
    case {1,3, 5}
        cfg = []; cfg.subjs_3 = S_data.subjs_3; cfg.subjs_2 = S_data.subjs_2;
        cfg.sFiles_3 = S_data.sFiles_3; cfg.sFiles_2 = S_data.sFiles_2;
        sub_demog_data = ecpfunc_read_sub_demog(cfg);
    case {2,4}
        clc
        cfg = []; cfg.subjs = S_data.subjs;
        cfg.sFiles = S_data.sFiles_32;
        sub_demog_data = ecpfunc_read_sub_demog_contrast(cfg);
end

%% TLE side (PT only)
patn_neuropsych_data = ecpfunc_read_patn_neuropsych();

%% Subject demog details
switch LI_analysis
    case {1,3 ,5}
        cfg = []; cfg.subjs_3 = S_data.subjs_3; cfg.subjs_2 = S_data.subjs_2;
        cfg.sFiles_3 = S_data.sFiles_3; cfg.sFiles_2 = S_data.sFiles_2;
        sub_demog_data = ecpfunc_read_sub_demog(cfg);
    case {2, 4, 6}
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
cfg = []; cfg.strt = 0; cfg.spt = 2; cfg.overlap = 0.01; cfg.linterval = 0.3;
wi  = do_time_intervals(cfg);

%% HCP atlas
glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat';
src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';

cfg = [];
cfg.src_fname = src_fname;
cfg.glass_dir = glass_dir;
cfg.glass_atlas = glass_atlas;
cfg.plotflag = 0;
Data_hcp_atlas = ecpfunc_hcp_atlas2(cfg);
net_sel_mutiple_label = Data_hcp_atlas.groups_labels';

%%
network_sel = [1:3,6:11];
colr = distinguishable_colors(length(network_sel));

%%
mLI_sub_hc = squeeze(nanmean(LI_hc.LI_sub,2));
mLI_sub_pt = squeeze(nanmean(LI_pt.LI_sub,2));

%%
% clc
patn_neuropsych_tle = ecpfunc_read_patn_neuropsych_tle();
TLESide = patn_neuropsych_tle.TLESide; SUBNO = patn_neuropsych_tle.SUBNO;

SUBNO_pt = [];
for i=1:length(LI_pt_ID)
    SUBNO_pt(i) = str2double(LI_pt_ID{i}(3:end));
end

[~,~,IB] = intersect(SUBNO_pt, SUBNO);
TLESide_sel = TLESide(IB);

TLE_left = find(TLESide_sel == 'Left');

LI_pt_val_left = LI_pt.LI_sub(:,TLE_left,:);

mLI_sub_left = squeeze(mean(LI_pt_val_left,2));

%%
LI_class_label = {'HC', 'PT'};

for j=1:length(LI_class_label)
    
    switch LI_class_label{j}
        case 'HC'
            d_in  = mLI_sub_hc;
        case 'PT'
            d_in  = mLI_sub_pt;
        case 'HC-PT'
            d_in  = mLI_sub_hc - mLI_sub_pt;
        case 'PT-Left'
            d_in  = mLI_sub_left;
        case 'HC - PT-Left'
            d_in  = mLI_sub_hc - mLI_sub_left;
    end
    
    cfg = [];
    cfg.data = d_in(network_sel, :);
    cfg.labels = net_sel_mutiple_label(network_sel);
    cfg.colors = colr;
    cfg.titleText = [LI_method_label{LI_method},', ', LI_class_label{j}];
    cfg.wi = wi;
    plotData(cfg);
    
    % - export figs
    cfg = []; cfg.outdir = save_dir; filename = LI_class_label{j}; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
end

%% Plot power - group mean
% close all
if  LI_method ==1
    
    pow_hc = transformPowSubTo3DArrays(LI_hc.pow_sub);
    pow_pt = transformPowSubTo3DArrays(LI_pt.pow_sub);
    
    mPow_sub_hc_left = 20.*squeeze(nanmean(pow_hc.left,2)); mPow_sub_hc_right = squeeze(nanmean(pow_hc.right,2));
    mPow_sub_pt_left = 20.*squeeze(nanmean(pow_pt.left,2)); mPow_sub_pt_right = squeeze(nanmean(pow_pt.right,2));
    
    cfg_pow = [];
    cfg_pow.labels = net_sel_mutiple_label(network_sel);
    cfg_pow.colors = colr;
    cfg_pow.wi = wi;
    cfg_pow.power_left = mPow_sub_hc_left(network_sel,:);
    cfg_pow.power_right = mPow_sub_hc_right(network_sel,:);
    cfg_pow.title = {'HC-left h', 'HC-right h'};
    plotPower(cfg_pow)
    
    % - export figs
    cfg = []; cfg.outdir = save_dir; filename = 'Power_hc'; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    
    cfg_pow.power_left = mPow_sub_pt_left(network_sel,:);
    cfg_pow.power_right = mPow_sub_pt_right(network_sel,:);
    cfg_pow.title = {'pt-left h', 'pt-right h'};
    plotPower(cfg_pow)
    
    % - export figs
    cfg = []; cfg.outdir = save_dir; filename = 'Power_pt'; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
end


%%
network_sel = 11;
network_sel = [1:3,6:11];
network_sel = [1,2,6,11];

% Assuming LI_method == 1 is still your condition to proceed
if LI_method == 1
    
    pow_hc = transformPowSubTo3DArrays(LI_hc.pow_sub);
    pow_pt = transformPowSubTo3DArrays(LI_pt.pow_sub);
    
    % Loop through each subject instead of calculating mean
    for subj = 1:1 % size(pow_hc.left, 2)
        cfg_pow = [];
        cfg_pow.labels = net_sel_mutiple_label(network_sel);
        cfg_pow.colors = colr;
        cfg_pow.wi = wi;
        
        % Use individual subject's power data for plotting
        cfg_pow.power_left = (pow_hc.left(network_sel,subj,:));
        cfg_pow.power_right = (pow_hc.right(network_sel,subj,:));
        cfg_pow.title = {['HC-left h - Subj ' num2str(subj)], ['HC-right h - Subj ' num2str(subj)]};
        plotPower(cfg_pow);
        
        % Export figs for each subject
        cfg = []; cfg.outdir = save_dir; filename = ['Power_hc_Subj_' num2str(subj)]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg);
    end
    
    % Loop through each subject instead of calculating mean
    for subj = 1:1 %size(pow_pt.left, 2)
        cfg_pow = [];
        cfg_pow.labels = net_sel_mutiple_label(network_sel);
        cfg_pow.colors = colr;
        cfg_pow.wi = wi;
        
        % Use individual subject's power data for plotting
        cfg_pow.power_left = (pow_pt.left(network_sel,subj,:));
        cfg_pow.power_right = (pow_pt.right(network_sel,subj,:));
        cfg_pow.title = {['HC-left h - Subj ' num2str(subj)], ['PT-right h - Subj ' num2str(subj)]};
        plotPower(cfg_pow);
        
        % Export figs for each subject
        cfg = []; cfg.outdir = save_dir; filename = ['Power_pt_Subj_' num2str(subj)]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg);
    end
end


