
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
LI_analysis_label = {'DICS_baseline','DICS_contrast','LCMV_basline','LCMV_contrast','DICS_anim'};

disp('1) DICS_baseline (Anim-vs-bsl vs. Symb-vs-bsl)')
disp('2) DICS_contrast (Anim-vs-Symb)')
disp('3) LCMV_basline (Anim-vs-bsl vs. Symb-vs-bsl)')
disp('4) LCMV_contrast (Anim-vs-Symb)')
disp('5) DICS_anim (Anim)')

LI_analysis = input('');

%%
LI_method_label = {'Magnitude', 'Counting','Bootstrapping'};

disp('1: Magnitude')
disp('2: Counting')
disp('3: Bootstrapping')
% LI_method = input(':');

%%
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/';
cd(data_save_dir)

%%
save_dir = fullfile(data_save_dir,LI_analysis_label{LI_analysis}, 'compare_LIs_Left_TLE');
checkOrCreateDir(save_dir)
cd(save_dir)

LI_pt = [];
for i=1:length(LI_method_label)
    LI_pt.(LI_method_label{i}) = load(fullfile(data_save_dir,LI_analysis_label{LI_analysis},LI_method_label{i},'LI_Patn'));
end

LI_pt.LI_methods = LI_method_label;

%% Run Brainstorm
Run_BS

%%
Run_load_surface_template

%% Read fMRI lats
fmri_LIs = ecpfunc_read_fmri_lat();

%%
datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';
BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_all_subjects';

cfg = [];
cfg.protocol = protocol;
cfg.datadir = datadir;
cfg.BS_data_dir = BS_data_dir;

switch LI_analysis
    case {1,5}
        cfg.datatag = 'wDICS_22_4_baseline';
        S_data = ecpfunc_read_sourcemaps_dics(cfg);
    case 2
        cfg.datatag = 'wDICS_contrast_18_4';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 3
        S_data = ecpfunc_read_sourcemaps(cfg);
    case 4
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

%%
switch LI_analysis
    case {1,3, 5}
        cfg = [];
        cfg.sub_demog_data = sub_demog_data;
        cfg.patn_neuropsych_data = patn_neuropsych_data;
        sub_TLE_sub_data = ecpfunc_read_sub_TLE_sub(cfg);
    case {2,4}
        clc
        cfg = [];
        cfg.sub_demog_data = sub_demog_data;
        cfg.patn_neuropsych_data = patn_neuropsych_data;
        sub_TLE_sub_data = ecpfunc_read_sub_TLE_sub_contrast(cfg);
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
    case {2,4}
        disp('1: Ctrl')
        disp('2: Patn')
        cfg = [];
        cfg.sub_demog_data = sub_demog_data;
        cfg.select_data = 1; S_data_hc = ecpfunc_select_contrast_data(cfg);
        cfg.select_data = 2; S_data_pt = ecpfunc_select_contrast_data(cfg);
        LI_pt_ID = S_data_pt.sFiles_subid;
end

LI_pt.pt_ID = LI_pt_ID;

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
% mLI_sub_hc = squeeze(nanmean(LI_hc.LI_sub,2));
% mLI_sub_pt = squeeze(nanmean(LI_pt.LI_sub,2));

%% L-TLE subjects
patn_neuropsych_tle = ecpfunc_read_patn_neuropsych_tle();

cfg = [];
cfg.patn_neuropsych_tle = patn_neuropsych_tle;
cfg.LI_pt = LI_pt;
LI_pt_left_tle = ecpfunc_select_left_tle(cfg);

%% MEG vs. fMRI lat analysis (PT)
[sub_MF_pt,IA,IB] = intersect(LI_pt_left_tle.pt_ID, fmri_LIs.ID.language_Lateral');

% LI_pt_val = LI_pt.(LI_method_label{1}).LI_sub;

% LI_pt_val_new = LI_pt_val(:,IA,:);
fmri_LIs_val = (fmri_LIs.val.language_Lateral(IB)); % Lateral regions was used.

% missing fMRI data.
difference = setdiff(LI_pt_left_tle.pt_ID, sub_MF_pt');
disp('missing from fMRI')
disp(difference');

%% MEG LI vs fMRI LI (language_Lateral)
% pause, close all,

LI_pt_val_new = [];
for i=1:length(LI_method_label)
    LI_pt_val_new.(LI_method_label{i}) = LI_pt_left_tle.(LI_method_label{i}).LI_sub(:, IA,:);
end

% Corr, MEG-fMRI
cfg = []; cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.ternary = 0;
cfg.thre = .2;
cfg.savefig = 1;
cfg.bf = 5;
cfg.outdir = save_dir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new;
cfg.fmri_LIs_val = fmri_LIs_val;
cfg.LI_method_label = LI_method_label;
% cfg.net_sel = [1,2,6];
cfg.net_sel = [11]; % 6, 1, 2
[megLI_sub_pt, fmri_LIs_val, ~] = do_MEG_fMRI_corr_contrast_all(cfg);

%% MEG LI vs fMRI LI (Ternary language_Lateral)
% pause, close all,
% clc

cfg = [];
cfg.thre = .2; cfg.LI = fmri_LIs_val;
fmri_LIs_trn = do_ternary_classification(cfg);
size(fmri_LIs_trn);

% concordance, MEG-fMRI
cfg = [];
cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.ternary = 1;
cfg.savefig = 1;
cfg.outdir = save_dir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new;
cfg.fmri_LIs_val = fmri_LIs_trn;
cfg.LI_method_label = LI_method_label;
% cfg.net_sel = [1,2,6];
cfg.net_sel = [11];
cfg.thre = 0.1;
cfg.buffervalue = 10;
[megLIs_trn, fmri_LIs_trn] = do_MEG_fMRI_concordance_contrast_all(cfg);

%% mean MEG li vs. fMRI
% pause, close all,

for i=1:length(LI_method_label)
    
    cfg = [];
    cfg.DataArray = [megLIs_trn{i}, fmri_LIs_trn];
    cfg.savefig = 1;
    cfg.outdir = save_dir;
    cfg.title = ['trn-', LI_method_label{i}];%'SD task, 58 PTs, raw LIs';
    do_plot_LIs(cfg)
    
    cfg = [];
    cfg.DataArray = [megLI_sub_pt{i}, fmri_LIs_val];
    cfg.savefig = 1;
    cfg.outdir = save_dir;
    cfg.title = ['LI-', LI_method_label{i}];%'SD task, 58 PTs, raw LIs';
    do_plot_LIs(cfg)
    %
    %     size(fmri_LIs_val)
    
    %
    [C, order] = confusionmat(megLIs_trn{i}, fmri_LIs_trn);
    % Visualize the confusion matrix
    figure;
    h = heatmap(C);
    h.XDisplayLabels = {'-1','0', '1'};
    h.YDisplayLabels = {'-1','0', '1'};
    xlabel('MEG');
    ylabel('fMRI');
    title(LI_method_label{i})%'Confusion Matrix');
    colorbar off
    colorbar('location', 'eastoutside')
    
    %%
    % - export figs
    cfg = []; cfg.outdir = save_dir; cfg.filename = ['Confusion Matrix: ', LI_method_label{i}];
    cfg.type = 'fig'; do_export_fig(cfg)

    
    % clc
    cfg = []; cfg.wi = wi;
    cfg.ID = sub_MF_pt;
    cfg.thre = 0.1;
    cfg.bf = 10;
    cfg.ternary = 0;
    cfg.savefig = 0;
    cfg.outdir = save_dir;
    cfg.net_sel_mutiple_label = net_sel_mutiple_label;
    cfg.net_sel_id = [1,2,5,6,11];
    cfg.lang_id = {'language_Angular'; 'language_Frontal';'language_PCingPrecun'; 'language_Temporal'; 'language_Lateral'};
    cfg.LI_val = LI_pt_val_new.(LI_method_label{i});
    cfg.fmri_LIs_val = fmri_LIs;
    cfg.idx = IB;
    cfg.title = LI_method_label{i};
    crr = do_MEG_fMRI_corr_contrast_rois(cfg);
    
    % - export figs
    cfg = []; cfg.outdir = save_dir; filename = ['net ROIs_', LI_method_label{i}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
end

cd(save_dir)