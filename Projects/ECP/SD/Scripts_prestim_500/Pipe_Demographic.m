%% ECP Semantic Decision Task Dataset, Medical College of Wisconsin

% Script: BS Process (Laterality Analysis)
% Project: ECP_SD
% Written by: Vahab Youssof Zadeh

clear; clc; close all; warning off;

%% Paths
restoredefaultpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/run')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/function')
Run_setpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other/')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/tools/helpful_tools/daviolinplot/daboxplot');

%% Define thresholds for MEG and fMRI
MEG_thre = 0.2; % MEG threshold
fMRI_thre = 0.1; % fMRI threshold

% Define the time interval bounds
lowerBound = 0.5;
upperBound = 1;

%% Laterality Analysis Labels
LI_analysis_label = {'DICS_indirect','DICS_directcontrast','LCMV_anim_vs_Symb','-','DICS_anim', 'DICS_contrast_prestim', 'dSPM_contrast'};
for i = 1:length(LI_analysis_label)
    disp([num2str(i) ') ' LI_analysis_label{i}]);
end
LI_analysis = input('');

%% Laterality Method Labels
LI_method_label = {'Magnitude', 'Counting','Bootstrapping', 'mean_combined'};
disp('1: Magnitude')
disp('2: Counting')
disp('3: Bootstrapping')
disp('4: Combined')
LI_method = input(':');

%% Data Save Directory
data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim';
cd(data_save_dir)

%% Processing Based on LI Method
switch LI_method
    case {1,2,3}
        save_dir = fullfile(data_save_dir, LI_analysis_label{LI_analysis}, LI_method_label{LI_method});
        checkOrCreateDir(save_dir)
        cd(save_dir)
        LI_hc = load('LI_Ctrl');
        LI_pt = load('LI_Patn');
    case 4
        for i = 1:3
            save_dir = fullfile(data_save_dir, LI_analysis_label{LI_analysis}, LI_method_label{i});
            checkOrCreateDir(save_dir)
            cd(save_dir)
            LI_hc = load('LI_Ctrl');
            LI_pt = load('LI_Patn');
            LI_hc_combined.(LI_method_label{i}) = LI_hc;
            LI_pt_combined.(LI_method_label{i}) = LI_pt;
        end
        save_dir = fullfile(data_save_dir, LI_analysis_label{LI_analysis}, LI_method_label{4});
        checkOrCreateDir(save_dir)
        cd(save_dir)
        LI_pt_combined.mean = (LI_pt_combined.(LI_method_label{1}).LI_sub + LI_pt_combined.(LI_method_label{2}).LI_sub + LI_pt_combined.(LI_method_label{3}).LI_sub) / 3;
        LI_hc_combined.mean = (LI_hc_combined.(LI_method_label{1}).LI_sub + LI_hc_combined.(LI_method_label{2}).LI_sub + LI_hc_combined.(LI_method_label{3}).LI_sub) / 3;
        LI_pt.LI_sub = LI_pt_combined.mean;
        LI_hc.LI_sub = LI_hc_combined.mean;
end

%% Brainstorm
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3';
BS_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/';
protocol = fullfile(BS_dir, 'data_full/protocol.mat');

Run_load_surface_template

%% fMRI LIs
fmri_LIs = ecpfunc_read_fmri_lat();

%% Directories and Configuration
datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';
BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_full';
cfg = [];
cfg.protocol = protocol;
cfg.datadir = datadir;
cfg.BS_data_dir = BS_data_dir;

%% Load Source Maps Based on Analysis Type
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

%% Subject Demographic Details
switch LI_analysis
    case {1,3,5}
        cfg = [];
        cfg.subjs_3 = S_data.subjs_3;
        cfg.subjs_2 = S_data.subjs_2;
        cfg.sFiles_3 = S_data.sFiles_3;
        cfg.sFiles_2 = S_data.sFiles_2;
        sub_demog_data = ecpfunc_read_sub_demog(cfg);
    case {2,4}
        clc
        cfg = [];
        cfg.subjs = S_data.subjs;
        cfg.sFiles = S_data.sFiles_32;
        sub_demog_data = ecpfunc_read_sub_demog_contrast(cfg);
end

%% TLE Side (PT Only)
patn_neuropsych_data = ecpfunc_read_patn_neuropsych();

%% Further Analysis Based on TLE Side
switch LI_analysis
    case {1,3,5}
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

%% Inter-Subject (Group) Averaging
switch LI_analysis
    case {1,3,5}
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

%% Time Intervals
cfg = []; cfg.strt = -0.5; cfg.spt = 2; cfg.overlap = 0.05; cfg.linterval = 0.3;
wi  = do_time_intervals(cfg);

%% HCP Atlas
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

%% Network Selection and Colors
network_sel = [1:3,6:11];
colr = distinguishable_colors(length(network_sel));

%% Mean Lateralization Index Calculation
mLI_sub_hc = squeeze(nanmean(LI_hc.LI_sub,2));
mLI_sub_pt = squeeze(nanmean(LI_pt.LI_sub,2));

%% TLE Side
patn_neuropsych_tle = ecpfunc_read_patn_neuropsych_tle();
TLESide = patn_neuropsych_tle.TLESide; SUBNO = patn_neuropsych_tle.SUBNO;
SUBNO_pt = [];
for i=1:length(LI_pt_ID)
    SUBNO_pt(i) = str2double(LI_pt_ID{i}(3:end));
end
[~,~,IB_tle] = intersect(SUBNO_pt, SUBNO);
TLESide_sel = TLESide(IB_tle);
TLE_left = find(TLESide_sel == 'Left');
LI_pt_val_left = LI_pt.LI_sub(:,TLE_left,:);
mLI_sub_left = squeeze(mean(LI_pt_val_left,2));

%% Demographic Details
patn_MEGfMRI_neuropsych = patn_neuropsych_tle(IB_tle,:);
summaryStats = @(x) table(mean(x,'omitnan'), median(x,'omitnan'), std(x,'omitnan'));
stats = varfun(summaryStats, patn_MEGfMRI_neuropsych, 'InputVariables', ...
    {'Age', 'EHQ', 'EducationYRS', 'WASI_BlckR_ZScore', 'WASI_VocR_ZScore'}, ...
    'OutputFormat', 'table');
disp(stats);

variablesToCorrelate = {'Age', 'EHQ', 'EducationYRS', 'WASI_BlckR_ZScore', 'WASI_VocR_ZScore'};
correlationMatrix = corr(patn_MEGfMRI_neuropsych{:, variablesToCorrelate}, 'Rows', 'complete');
disp(array2table(correlationMatrix, 'VariableNames', variablesToCorrelate, 'RowNames', variablesToCorrelate));

lm = fitlm(patn_MEGfMRI_neuropsych, 'Age~EHQ');
disp(lm);

% Gender summary
genderSummary = groupsummary(patn_MEGfMRI_neuropsych, 'Sex', 'IncludeMissingGroups', true);
disp(genderSummary);

% TLE side summary
TLESummary = groupsummary(patn_MEGfMRI_neuropsych, 'TLESide', 'IncludeMissingGroups', true);
disp(TLESummary);

% Age statistics
meanAge = mean(patn_MEGfMRI_neuropsych.Age, 'omitnan');
stdAge = std(patn_MEGfMRI_neuropsych.Age, 'omitnan');
disp(['Mean Age: ', num2str(meanAge)]);
disp(['Standard Deviation of Age: ', num2str(stdAge)]);

minAge = min(patn_MEGfMRI_neuropsych.Age);
maxAge = max(patn_MEGfMRI_neuropsych.Age);
ageRange = maxAge - minAge;
disp(['Minimum Age: ', num2str(minAge)]);
disp(['Maximum Age: ', num2str(maxAge)]);
disp(['Age Range: ', num2str(ageRange)]);

% Handedness summary
handednessSummary = groupsummary(patn_MEGfMRI_neuropsych, 'DomntHand', 'IncludeMissingGroups', true);
disp(handednessSummary);
