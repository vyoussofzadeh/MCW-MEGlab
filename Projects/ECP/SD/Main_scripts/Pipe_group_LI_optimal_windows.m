
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
% LI_method = input(':');

%%
data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim';
cd(data_save_dir)

%%
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
    case {1,5}
        cfg.datatag = 'wDICS_baseline_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics(cfg);
    case 2
        cfg.datatag = 'wDICS_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 3
        cfg.datamask = fullfile('./Group_analysis/LCMV/results_average*.mat');
        S_data = ecpfunc_read_sourcemaps(cfg);
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
    case {2,4,6}
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
    case {2,4,6}
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
cfg = []; cfg.strt = -0.5; cfg.spt = 2; cfg.overlap = 0.05; cfg.linterval = 0.3;
wi  = do_time_intervals(cfg);

%% HCP atlas
glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat';
% glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_HCP_MMP1_360_2024.mat';
src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';

cfg = [];
cfg.src_fname = src_fname;
cfg.glass_dir = glass_dir;
cfg.glass_atlas = glass_atlas;
cfg.plotflag = 1;
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
            rSNR_new.Bootstrapping(i, j).left = mean(rSNR_new.Bootstrapping(i, j).left, 2);
        end
        
        % Calculate mean for Bootstrapping right if necessary
        if size(rSNR_new.Bootstrapping(i, j).right, 2) > 1
            rSNR_new.Bootstrapping(i, j).right = mean(rSNR_new.Bootstrapping(i, j).right, 2);
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
% clc

plot_indiv_LI = 0; plot_rSNR = 0; plot_rSNR_LI = 0;

% lowerBound = 0.3; upperBound = 1.3;
% lowerBound = 0.45; upperBound = 1.2;
lowerBound = 0.5; upperBound = 1.5;
lowerBound = 0.5; upperBound = 1.6;
lowerBound = 0.90; upperBound = 1.6;

% lowerBound = 0.45; upperBound = 0.75;


% opt_method = 'LI'; %'rsnr'
opt_method = 'rsnr';
% opt_method = 'AUC';
% opt_method = 'DomH';
% opt_method = 'rsnrmax';
% opt_method = 'wrsnr';

run_optimalLIs_snr_dics
% run_optimalLIs

%%
cfg  = [];
cfg.fmri_LIs_val = fmri_LIs_val; cfg.wi = wi;
cfg.sub_MF_pt = sub_MF_pt;
cfg.LI_method_label = LI_method_label;
cfg.LI_pt_val_new = LI_pt_val_new;
cfg.fmri_LIs = fmri_LIs;
cfg.IB_megfmri = IB_megfmri;
cfg.rSNR_new = rSNR_new;
cfg.opt_method = opt_method; cfg.lowerBound = lowerBound; cfg.upperBound = upperBound;
cfg.MEG_thre = MEG_thre; cfg.fMRI_thre = fMRI_thre;
cfg.LI_method_labels = LI_method_labels;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.plot_indiv_LI = plot_indiv_LI;
cfg.plot_rSNR = plot_rSNR; cfg.plot_rSNR_LI = plot_rSNR_LI;
cfg.idcx =  [1, 2, 6, 11];
MConcordance = func_run_optimalLIs_snr_dics(cfg);

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

%% Response (Reaction) Time Data
load('/data/MEG/Research/aizadi/process/RT_summary/ResponseTime.mat')
[~,~,IB_reactiontime] = intersect(sub_MF_pt, T.Sub_ID);
T_patn_MEGfMRI = T(IB_reactiontime,:);
meanAnimal = mean(T_patn_MEGfMRI.Animal, 'omitnan');
stdAnimal = std(T_patn_MEGfMRI.Animal, 'omitnan');
meanSymbol = mean(T_patn_MEGfMRI.Symbol, 'omitnan');
stdSymbol = std(T_patn_MEGfMRI.Symbol, 'omitnan');
disp(['Mean of Animal reaction times: ', num2str(meanAnimal)]);
disp(['Standard Deviation of Animal reaction times: ', num2str(stdAnimal)]);
disp(['Mean of Symbol reaction times: ', num2str(meanSymbol)]);
disp(['Standard Deviation of Symbol reaction times: ', num2str(stdSymbol)]);
[correlationCoefficient, p] = corr(T_patn_MEGfMRI.Animal, T_patn_MEGfMRI.Symbol, 'Rows', 'complete');
validPairs = sum(~isnan(T_patn_MEGfMRI.Animal) & ~isnan(T_patn_MEGfMRI.Symbol));
df = validPairs - 2;
disp(['Correlation coefficient: ', num2str(correlationCoefficient)]);
disp(['Degrees of freedom: ', num2str(df)]);

%% Task Performance
taskperf_datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/';
sub_MF_pt_num = cellfun(@(x) str2double(x(3:end)), sub_MF_pt);
taskPerformanceDataPath = fullfile(taskperf_datadir, 'TaskPerformanceSD.mat');
load(taskPerformanceDataPath); 
[~,~,IB_taskperformance] = intersect(sub_MF_pt_num, accuracyResults.Subject);
accuracyResults_updt = accuracyResults(IB_taskperformance,:);
meanAccBySubject_Animal = groupsummary(accuracyResults_updt, 'Subject', 'mean', 'Animal_ACC');
meanAccBySubject_Falsefont = groupsummary(accuracyResults_updt, 'Subject', 'mean', 'Falsefont_ACC');
totalmean = mean(meanAccBySubject_Falsefont.mean_Falsefont_ACC);
disp(['The total mean of mean_Falsefont_ACC is: ', num2str(totalmean)]);
totalstd = std(meanAccBySubject_Falsefont.mean_Falsefont_ACC);
disp(['The total std of mean_Falsefont_ACC is: ', num2str(totalstd)]);
totalmean = mean(meanAccBySubject_Animal.mean_Animal_ACC);
disp(['The total mean of mean_Anim_ACC is: ', num2str(totalmean)]);
totalstd = std(meanAccBySubject_Animal.mean_Animal_ACC);
disp(['The total std of mean_Anim_ACC is: ', num2str(totalstd)]);

%% Investigate Discordant Samples of Best Results and Obtain Corresponding MEG_LI and fMRI_LI
run_plot_bestLIs
