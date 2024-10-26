
%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 05/21/2024

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/run')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/function')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/tools/helpful_tools/daviolinplot/daboxplot')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Main_scripts/run')
Run_setpath

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other/')

%%
MEG_thre = 10; % MEG threshold
fMRI_thre = 10; % fMRI threshold

figtype = 'svg';

%%
LI_analysis_label = {'DICS_indirect','DICS_directcontrast','LCMV_anim_vs_Symb','-','DICS_anim', 'DICS_contrast_prestim', 'dSPM_contrast','DICS_directcontrast_localthresh'};

for i = 1:length(LI_analysis_label)
    disp([num2str(i) ') ' LI_analysis_label{i}]);
end

LI_analysis = input('');

%%
LI_method_label = {'Magnitude', 'Counting','Bootstrapping'};

disp('1: Magnitude')
disp('2: Counting')
disp('3: Bootstrapping')
% LI_method = input(':');

%%
data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim';
cd(data_save_dir)

%%
clc
save_dir = fullfile(data_save_dir,LI_analysis_label{LI_analysis}, 'compare_LIs');
checkOrCreateDir(save_dir)
cd(save_dir)

LI_pt = []; rSNR = [];

for i=1:length(LI_method_label)
    LI_pt.(LI_method_label{i}) = load(fullfile(data_save_dir,LI_analysis_label{LI_analysis},LI_method_label{i},'LI_Patn'));
    if i == 1
        rSNR.(LI_method_label{i}) = LI_pt.(LI_method_label{i}).pow_sub;
    else
        rSNR.(LI_method_label{i}) = LI_pt.(LI_method_label{i}).count_sub;
    end
end

%% Run Brainstorm
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3';

BS_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/';
protocol = fullfile(BS_dir, 'data_full/protocol.mat');

%%
Run_load_surface_template

%% Read fMRI lats
fmri_LIs = ecpfunc_read_fmri_lat();

%%
datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';
BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_full';

cfg = [];
cfg.protocol = protocol;
cfg.datadir = datadir;
cfg.BS_data_dir = BS_data_dir;

switch LI_analysis
    case 1
        cfg.datatag = 'wDICS_baseline_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics(cfg);
    case {2,8}
        cfg.datatag = 'wDICS_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 3
        cfg.datamask = fullfile('./Group_analysis/LCMV/results_average*.mat');
        S_data = ecpfunc_read_sourcemaps(cfg);
    case 5
        protocol = fullfile(BS_dir, 'data/protocol.mat');
        BS_data_dir = fullfile(BS_dir,'data');
        cfg.protocol = protocol;
        cfg.BS_data_dir = BS_data_dir;
        cfg.datatag = 'wDICS_baseline_18_4';
        S_data = ecpfunc_read_sourcemaps_dics_anim(cfg);
    case 6
        cfg.datatag = 'wDICS_contrast_18_4_50ms';
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
    case {2,4,6,8}
        %         clc
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
    case {2,4,6,8}
        %         clc
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
    case {2,4,8}
        disp('1: Ctrl')
        disp('2: Patn')
        cfg = [];
        cfg.sub_demog_data = sub_demog_data;
        cfg.select_data = 1; S_data_hc = ecpfunc_select_contrast_data(cfg);
        cfg.select_data = 2; S_data_pt = ecpfunc_select_contrast_data(cfg);
        LI_pt_ID = S_data_pt.sFiles_subid;
end

%%
cfg = []; cfg.strt = -0.5; cfg.spt = 2; cfg.overlap = 0.05; cfg.linterval = 0.3;
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

net_sel_mutiple_label{1} = 'Ang';
net_sel_mutiple_label{2} = 'Front';
net_sel_mutiple_label{6} = 'Temp';
net_sel_mutiple_label{11} = 'Lat';

%%
network_sel = [1:3,6:11];
colr = distinguishable_colors(length(network_sel));

%%
patn_neuropsych_tle = ecpfunc_read_patn_neuropsych_tle();
TLESide = patn_neuropsych_tle.TLESide; SUBNO = patn_neuropsych_tle.SUBNO;

SUBNO_pt = [];
for i=1:length(LI_pt_ID)
    SUBNO_pt(i) = str2double(LI_pt_ID{i}(3:end));
end

[~,~,IB] = intersect(SUBNO_pt, SUBNO);
TLESide_sel = TLESide(IB);

TLE_left = find(TLESide_sel == 'Left');

%% MEG vs. fMRI lat analysis (PT)
[sub_MF_pt,IA,IB_megfmri] = intersect(LI_pt_ID, fmri_LIs.ID.language_Lateral');

LI_pt_val = LI_pt.(LI_method_label{1}).LI_sub;

fmri_LIs_val = (fmri_LIs.val.language_Lateral(IB_megfmri)); % Lateral regions was used.

% missing from fMRI
difference = setdiff(LI_pt_ID, sub_MF_pt');
disp('missing from fMRI')
disp(difference');

%% MEG LI vs fMRI LI (language_Lateral)
LI_pt_val_new = [];
rSNR_new = [];

for i=1:length(LI_method_label)
    LI_pt_val_new.(LI_method_label{i}) = LI_pt.(LI_method_label{i}).LI_sub(:, IA,:);
    rSNR_new.(LI_method_label{i}) = rSNR.(LI_method_label{i})(:, IA);
end

cfg = [];
cfg.thre = fMRI_thre; cfg.LI = fmri_LIs_val;
fmri_LIs_trn = do_ternary_classification2(cfg);
size(fmri_LIs_trn);

Pt_ID  = LI_pt.Magnitude.setup.S_data.sFiles_subid(IA)';
Pt_sfiles  = LI_pt.Magnitude.setup.S_data.sFiles_in(IA)';


%%
% Loop through each element in the struct arrays
for i = 1:size(rSNR_new.Magnitude, 1)
    for j = 1:size(rSNR_new.Magnitude, 2)
        % Adjust Counting to make it a column vector
        if isrow(rSNR_new.Counting(i, j).left)
            rSNR_new.Counting(i, j).left = rSNR_new.Counting(i, j).left.';
        end
        
        % Adjust Counting for right if necessary
        if isrow(rSNR_new.Counting(i, j).right)
            rSNR_new.Counting(i, j).right = rSNR_new.Counting(i, j).right.';
        end

        % Calculate mean for Bootstrapping if it's not a single column vector
        if size(rSNR_new.Bootstrapping(i, j).left, 2) > 1
            rSNR_new.Bootstrapping(i, j).left = nanmean(rSNR_new.Bootstrapping(i, j).left, 2);
        end
        
        % Calculate mean for Bootstrapping right if necessary
        if size(rSNR_new.Bootstrapping(i, j).right, 2) > 1
            rSNR_new.Bootstrapping(i, j).right = nanmean(rSNR_new.Bootstrapping(i, j).right, 2);
        end

        % Ensure Magnitude is a column vector for both left and right
        if isrow(rSNR_new.Magnitude(i, j).left)
            rSNR_new.Magnitude(i, j).left = rSNR_new.Magnitude(i, j).left.';
        end
        if isrow(rSNR_new.Magnitude(i, j).right)
            rSNR_new.Magnitude(i, j).right = rSNR_new.Magnitude(i, j).right.';
        end
    end
end

%% Plot MEG LI for selected network ROIs
run_plot_MEGLIs

run_plot_MEGLIs_ver2

%%
run_plot_MEGsnr

%% Fixed interval analysis
%- Summerize LIs_fixedinterva
run_sumLIs_fixedinterval

%- plot_fixedinterval_corcon_LIsmethods
run_plot_fixedinterval_corcon_LIsmethods

%- plot_fixedinterval_corcon_rois1
run_plot_fixedinterval_corcon_rois1

%- plot_fixedinterval_corcon_rois2
% run_plot_fixedinterval_corcon_rois2

%- plot_fixedinterval_rois
run_plot_fixedinterval_rois

%- table_fixedinterval
run_table_fixedinterval

%- plot_fixedinterval
run_plot_fixedinterval

%% Dynamic interval analysis
plot_indiv_LI = 0; plot_rSNR = 0; plot_rSNR_LI = 0;

disp('1) res-snr and optimal bound selection')
disp('2) res-snr with fixed lower and upper bounds')
disp('3) res-snr with optimzation against fMRI LIs - only for inspection!')
opt_sel = input('optimal method: ');

switch opt_sel
    case 1
        minlowerband = 0.35; maxUpperband = 1.5;
        opt_method = 'rsnr_optbound'; %'rsnr_optbound_mean'; 
        run_optimalLIs_snr_dics_rois
        disp(['mean_conc:', num2str(mean(summaryTableDynamic_save.Concordance))])
    case 2
        lowerBound = 0.35; upperBound = 1.5;
        opt_method = 'rsnr'; run_optimalLIs_snr_dics
    case 3
        opt_method = 'rsnr'; run_optimalLIs_usingfMRI
end

%%
%- plot_compare_fixed_opt_methods
run_plot_compare_fixed_opt_methods

%- plot_compare_fixed_opt_rois
run_plot_compare_fixed_opt_rois

%- Plot_compare_fixed_opt
run_plot_compare_fixed_opt

%- plot_compare_fixed_opt_rois_summ
run_plot_compare_fixed_opt_rois_summ
run_plot_compare_fixed_opt_gsum

%% Pick the Best Results Out of 3 LI Methods for Discordant Analyses
run_table_bestLIs

%% Response (reaction) time
run_responsereaction

%% Task Performance
run_taskperformance

%% Investigate Discordant Samples of Best Results and Obtain Corresponding MEG_LI and fMRI_LI
run_plot_bestLIs
