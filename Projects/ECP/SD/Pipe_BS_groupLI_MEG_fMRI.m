
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
LI_analysis_label = {'DICS_baseline','DICS_contrast','LCMV_basline','LCMV_contrast'};

disp('1) DICS_baseline (Anim-vs-bsl vs. Symb-vs-bsl)')
disp('2) DICS_contrast (Anim-vs-Symb)')
disp('3) LCMV_basline (Anim-vs-bsl vs. Symb-vs-bsl)')
disp('4) LCMV_contrast (Anim-vs-Symb)')

LI_analysis = input('');

%%
LI_method_label = {'Magnitude', 'Counting','Bootstrapping'};

disp('1: Magnitude')
disp('2: Counting')
disp('3: Bootstrapping')
LI_method = input(':');

%%
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/';
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
BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_all_subjects';

cfg = [];
cfg.protocol = protocol;
cfg.datadir = datadir;
cfg.BS_data_dir = BS_data_dir;

switch LI_analysis
    case 1
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
    case {1,3}
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
    case {1,3}
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
    case {1,3}
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

%%
cfg = []; cfg.strt = 0; cfg.spt = 2; cfg.overlap = 0.01; cfg.linterval = 0.3;
wi  = do_time_intervals(cfg);

%% Atlas
clc
glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat';
cfg = []; 
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
cfg.glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';
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
clc
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

% LI_pt_val_left = LI_pt.LI_sub(:,TLE_left,:);
% mLI_sub_left = squeeze(mean(LI_pt_val_left,2));

%%
clc
LI_class_label = {'HC', 'PT','HC-PT','PT-Left','HC - PT-Left'};

disp('1: HC')
disp('2: PT')
disp('3: HC-PT')
disp('4: PT-Left')
disp('5: HC - PT-Left')
% LI_class = input(':');


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
    plotData(cfg);
    
    % - export figs
    cfg = []; cfg.outdir = save_dir; filename = LI_class_label{j}; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
end

%% Power analysis
% switch LI_method
%     case 'Magnitude'
%         run_power_analysis
% end

%% MEG vs. fMRI lat analysis (PT)
clc
[sub_MF_pt,IA,IB] = intersect(LI_pt_ID, fmri_LIs.ID.language_Lateral');

LI_pt_val = LI_pt.LI_sub;

LI_pt_val_new = LI_pt_val(:,IA,:);
fmri_LIs_val = (fmri_LIs.val.language_Lateral(IB)); % Lateral regions was used.

% missing from fMRI
difference = setdiff(LI_pt_ID, sub_MF_pt');
disp('missing from fMRI')
disp(difference');

%% MEG LI vs fMRI LI (language_Lateral)
% close all
clc

% Corr, MEG-fMRI
cfg = []; cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.ternary = 0;
cfg.thre = .2;
cfg.savefig = 1;
cfg.bf = 10;
cfg.outdir = outdir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new; 
cfg.fmri_LIs_val = fmri_LIs_val; 
% cfg.net_sel = [1,2,6];
cfg.net_sel = [11]; % 6, 1, 2
[megLI_sub_pt, fmri_LIs_val, ~] = do_MEG_fMRI_corr_contrast(cfg);

%% MEG LI vs fMRI LI (Ternary language_Lateral)
% close all, 
% clc

cfg = [];
cfg.thre = .3; cfg.LI = fmri_LIs_val;
fmri_LIs_trn = do_ternary_classification(cfg);
size(fmri_LIs_trn);

% concordance, MEG-fMRI
cfg = [];
cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.ternary = 1;
cfg.savefig = 1;
cfg.outdir = outdir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new; 
cfg.fmri_LIs_val = fmri_LIs_trn; 
% cfg.net_sel = [1,2,6];
cfg.net_sel = [11];
cfg.thre = 0.15;
cfg.buffervalue = 4;
[megLIs_trn, fmri_LIs_trn] = do_MEG_fMRI_concordance_contrast(cfg);

% disp([megLIs_trn, fmri_LIs_trn])

idx = find(abs(megLIs_trn - fmri_LIs_trn) > 0);
sub_MF_pt(idx)

%%
% close all
% clc

cfg = [];
cfg.DataArray = [megLIs_trn, fmri_LIs_trn]; 
cfg.savefig = 1;
cfg.outdir = outdir; 
cfg.title = 'SD task, 58 PTs, ternary class';
do_plot_LIs(cfg)

cfg = [];
cfg.DataArray = [megLI_sub_pt, fmri_LIs_val]; 
cfg.savefig = 1;
cfg.outdir = outdir; 
cfg.title = 'SD task, 58 PTs, raw LIs';
do_plot_LIs(cfg)

size(fmri_LIs_val)

%%
[C, order] = confusionmat(megLIs_trn, fmri_LIs_trn);
% Visualize the confusion matrix
figure;
h = heatmap(C);
h.XDisplayLabels = {'-1','0', '1'};
h.YDisplayLabels = {'-1','0', '1'};
xlabel('MEG');
ylabel('fMRI');
title('Confusion Matrix');

% - export figs
cfg = []; cfg.outdir = outdir; cfg.filename = ['Confusion Matrix']; 
cfg.type = 'fig'; do_export_fig(cfg)

cd(outdir)

%% Corr. ROIs
fmri_LIs_val = (fmri_LIs.val.language_Frontal(IB)); net_sel = 2;
% fmri_LIs_val = (fmri_LIs.val.language_Temporal(IB)); net_sel = 6;
fmri_LIs_val = (fmri_LIs.val.language_Lateral(IB)); net_sel = 11;

cfg = []; cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.thre = 0.1;
cfg.bf = 10;
cfg.ternary = 0;
cfg.savefig = 0;
cfg.outdir = outdir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new;
size(LI_pt_val_new)
size(fmri_LIs_val)
cfg.fmri_LIs_val = fmri_LIs_val; cfg.net_sel = net_sel;
crr = do_MEG_fMRI_corr_contrast(cfg);

%% ROIs
clc
cfg = []; cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.thre = 0.1;
cfg.bf = 10;
cfg.ternary = 0;
cfg.savefig = 0;
cfg.outdir = outdir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new; 
cfg.fmri_LIs_val = fmri_LIs;
cfg.idx = IB;
crr = do_MEG_fMRI_corr_contrast_rois(cfg);

%% Corr.(tern)
% clc, close all

fmri_LIs_val = (fmri_LIs.val.language_Frontal(IB)); net_sel = 2;
fmri_LIs_val = (fmri_LIs.val.language_Angular(IB)); net_sel = 1;
fmri_LIs_val = (fmri_LIs.val.language_Temporal(IB)); net_sel = 6;
fmri_LIs_val = (fmri_LIs.val.language_Lateral(IB)); net_sel = 11;

cfg = [];
cfg.thre = .3; cfg.LI = fmri_LIs_val;
fmri_LIs_trn = do_ternary_classification(cfg);
size(fmri_LIs_trn);

%- concordance (similarity)
cfg = [];
cfg.wi = wi;
cfg.thre = 0.15;
cfg.ternary = 1;
cfg.savefig = 0;
cfg.outdir = outdir;
cfg.ID = sub_MF_pt;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new;
cfg.fmri_LIs_val = fmri_LIs_trn;
cfg.net_sel = net_sel;
cfg.buffervalue = 10;
conc = do_MEG_fMRI_concordance_contrast(cfg);

%%