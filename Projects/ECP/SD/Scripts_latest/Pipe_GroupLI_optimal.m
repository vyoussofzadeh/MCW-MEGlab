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
% MEG_thre = 0.2; % MEG threshold
% fMRI_thre = 0.1; % fMRI threshold

MEG_thre = 10; % MEG threshold
fMRI_thre = 10; % fMRI threshold

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

%% BS
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

%% Lateralization Index Analysis
LI_class_label = {'HC', 'PT', 'HC-PT', 'PT-Left', 'HC-PT-Left'};
for j = 1:2
    switch LI_class_label{j}
        case 'HC'
            d_in  = mLI_sub_hc;
        case 'PT'
            d_in  = mLI_sub_pt;
        case 'HC-PT'
            d_in  = mLI_sub_hc - mLI_sub_pt;
        case 'PT-Left'
            d_in  = mLI_sub_left;
        case 'HC-PT-Left'
            d_in  = mLI_sub_hc - mLI_sub_left;
    end
    cfg = [];
    cfg.data = d_in(network_sel, :);
    cfg.labels = net_sel_mutiple_label(network_sel);
    cfg.colors = colr;
    cfg.titleText = [LI_method_label{LI_method},', ', LI_class_label{j}];
    cfg.wi = wi;
    plotData(cfg);
    cfg = []; cfg.outdir = save_dir; filename = LI_class_label{j}; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
    combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
end

close all

%% MEG vs. fMRI Lat Analysis (PT)
[sub_MF_pt,IA,IB_megfmri] = intersect(LI_pt_ID, fmri_LIs.ID.language_Lateral');
LI_pt_val = LI_pt.LI_sub;
MEG_LI_Data = LI_pt_val(:,IA,:);
fMRI_LI = fmri_LIs.val.language_Lateral(IB_megfmri);

% Missing from fMRI
difference = setdiff(LI_pt_ID, sub_MF_pt');
disp('missing from fMRI')
disp(difference');

%- Demographic Details
patn_MEGfMRI_neuropsych = ecpfunc_sum_neuropsych(sub_MF_pt);

%% MEG LI vs fMRI LI (language_Lateral)
cfg = [];
cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.ternary = 1;
cfg.thre = 10;
cfg.savefig = 0;
cfg.bf = 1;
cfg.outdir = save_dir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = MEG_LI_Data;
cfg.fmri_LIs_val = fMRI_LI;
cfg.net_sel = [11];
[megLI_sub_pt, fmri_LIs_val, ~, interval_idx] = do_MEG_fMRI_corr_contrast2(cfg);

cfg = [];
cfg.thre = 10;
cfg.LI = fmri_LIs_val;
fmri_LIs_trn = do_ternary_classification2(cfg);
size(fmri_LIs_trn);

%% Constant Intervals LIs
MEG_LI = squeeze(MEG_LI_Data(11,:,:));
timePoints = mean(wi,2);
IntervalSize = 1;
[optimalInterval, correlations] = findOptimalMEGInterval(MEG_LI, fMRI_LI, timePoints, IntervalSize);
sub_IDs = sub_MF_pt;
nsub_IDs = cellfun(@(x) [num2str(find(strcmp(sub_IDs, x))), ':', x], sub_IDs, 'UniformOutput', false);
optimalInterval_constant = ones(size(MEG_LI, 1),1).*optimalInterval;
bf = 0.15;
lowerBound_constant = optimalInterval - bf;
upperBound_constant =  optimalInterval + bf;
[groupCorrelation_cnst, optimalTimePoints_cnst] = computeGroupLevelMEGfMRICorrelation_timepoints(MEG_LI, fMRI_LI, timePoints, lowerBound_constant, upperBound_constant);
[concordance_cnst, discordantSubs_cnst] = ...
    calculateConcordanceForTimePoints(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalInterval_constant);
disp(nsub_IDs(discordantSubs_cnst)')
subjectForPlot = discordantSubs_cnst;
findIndividualOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, subjectForPlot, lowerBound_constant, upperBound_constant);
cfg = []; cfg.outdir = save_dir; filename = ['ConstantTimePoints_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

close all

%% Optimal Time Points LIs
lowerBound = 0.4;
upperBound = 0.9;

close all
clc
[groupCorrelation, optimalTimePoints] = computeGroupLevelMEGfMRICorrelation_timepoints(MEG_LI, fMRI_LI, timePoints, lowerBound, upperBound);
% plotOptimalTimePointsOnMEG(MEG_LI, fMRI_LI, timePoints, optimalTimePoints, []);
meanOptimalTime = mean(optimalTimePoints);
disp(['Mean of optimal time points: ', num2str(meanOptimalTime)]);
[concordance, discordantSubs] = calculateConcordanceForTimePoints(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimePoints);
disp(['Concordance: ', num2str(concordance)]);
disp(['Correlation: ', num2str(groupCorrelation)]);

if ~isempty(discordantSubs)
    disp('Discordant Subjects:');
    disp(nsub_IDs(discordantSubs));
else
    disp('No Discordant Subjects Found');
end

%%
% plotOptimalTimePointsOnMEG(MEG_LI, fMRI_LI, [], timePoints, optimalTimePoints);
plotOptimalTimePointsOnMEG2(MEG_LI, fMRI_LI, timePoints, optimalTimePoints, discordantSubs, MEG_thre, lowerBound, upperBound);
cfg = []; cfg.outdir = save_dir; filename = ['Optimal_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

%% Mean Based Analysis
[groupCorrelation_mean] = computeGroupLevelMEGfMRICorrelation_avgLI(MEG_LI, fMRI_LI, timePoints, lowerBound, upperBound);
[concordance_mean, discordantSubs_mean] = calculateConcordanceUsingAvgMEGLI(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, lowerBound, upperBound);
disp(['Concordance_mean: ', num2str(concordance_mean)]);
disp(['Correlation_mean: ', num2str(groupCorrelation_mean)]);
findIndividualOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, discordantSubs_mean, lowerBound, upperBound);
cfg = []; cfg.outdir = save_dir; filename = ['mean_optimalTimePoints_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

close all

%% Power Values
if LI_method == 1
    
    % first method
    pow_hc = transformPowSubTo3DArrays(LI_hc.pow_sub);
    pow_pt = transformPowSubTo3DArrays(LI_pt.pow_sub);
    mPow_sub_hc_left = squeeze(nanmean(pow_hc.left,2)); mPow_sub_hc_right = squeeze(nanmean(pow_hc.right,2));
    mPow_sub_pt_left = squeeze(nanmean(pow_pt.left,2)); mPow_sub_pt_right = squeeze(nanmean(pow_pt.right,2));
    power_left = 20.*squeeze(pow_pt.left(11,IA,:));
    power_right = 20.*squeeze(pow_pt.right(11,IA,:));
    plot_flag = 0;
    
    [optimalTimePoints_pow] = plotSubjectPowerOverlay(power_left, power_right, sub_IDs, timePoints, plot_flag, lowerBound, upperBound);
    if plot_flag == 1
        cfg = []; cfg.outdir = save_dir; filename = ['Power_values_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
        close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    end
    [concordance_pow, discordantSubs_pow] = calculateConcordanceForTimePoints(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimePoints_pow);
    disp(['Concordance: ', num2str(concordance_pow)]);
    if ~isempty(discordantSubs_pow)
        disp('Discordant Subjects:');
        disp(nsub_IDs(discordantSubs_pow));
    else
        disp('No Discordant Subjects Found');
    end
    
    % Second method
    % - diff (L - R)
    [optimalTimePoints, interval] = plotSubjectPowerOverlay(power_left, power_right, sub_IDs, timePoints, 1, 0.2, 2.0);
    [optimalTimePoints, ~] = plotSubjectPowerOverlay(power_left, power_right, sub_IDs, timePoints, 0, interval(1), interval(2));
    
    % - Left vs. right
    [optimalTimePoints, interval] = plotSubjectPowerOverlay2(power_left, power_right, sub_IDs, timePoints, 1, 0.3, 1.5);
    [optimalTimePoints, ~] = plotSubjectPowerOverlay2(power_left, power_right, sub_IDs, timePoints, 0, interval(1), interval(2));
    
    % sum of left and right
    [optimalTimePoints, interval] = plotSubjectPowerOverlay3(power_left, power_right, sub_IDs, timePoints, 0, 0.4, 0.9);
    [optimalTimePoints, ~] = plotSubjectPowerOverlay2(power_left, power_right, sub_IDs, timePoints, 0, interval(1), interval(2));
    [groupCorrelation_cnst, optimalTimePoints] = computeGroupLevelMEGfMRICorrelation_timepoints(MEG_LI, fMRI_LI, timePoints, interval(1), interval(2));
    
    % - Area under the curve
    aucDiff = plotSubjectPowerOverlay4(power_left, power_right, sub_IDs, timePoints, 0,  0.1, 0.8);
    groupCorrelation = corr(aucDiff, fMRI_LI, 'Rows', 'complete')   
    
    if plot_flag == 1
        cfg = []; cfg.outdir = save_dir; filename = ['Power_values_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    end
    
    [concordance_pow, discordantSubs, MEG_LI_Sel] = ...
        calculateConcordanceForTimePoints(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimePoints);
    groupCorrelation = corr(MEG_LI_Sel', fMRI_LI, 'Rows', 'complete')
    
    %     [concordance_pow, discordantSubs_pow] = calculateConcordanceForTimePoints(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimePoints_pow);
    disp(['Concordance: ', num2str(concordance_pow)]);
    if ~isempty(discordantSubs)
        disp('Discordant Subjects:');
        disp(nsub_IDs(discordantSubs));
    else
        disp('No Discordant Subjects Found');
    end
    
    disp(nsub_IDs(discordantSubs)')
    subjectForPlot = discordantSubs;
%     [optimalTimePoints, LI_opt] = findIndividualOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, subjectForPlot, interval(1), interval(2));
    cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    if LI_method == 1
        plotIndividualSourcePower(power_left, power_right, MEG_LI, fMRI_LI, sub_IDs, timePoints, MEG_LI_Sel, subjectForPlot)
        cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_Dicordant_Pow_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    end    
end

cd(save_dir)

%% Optimal Time, Scatter Plot
numSubjects = length(optimalTimePoints);
subjectIDs = 1:numSubjects;
figure;
scatter(subjectIDs, optimalTimePoints, 'filled');
title('Opt Time');
xlabel('Subject ID');
ylabel('Optimal Time Point (s)');
grid on;
set(gca,'color','none');
set(gcf, 'Position', [400   100   300   300]);
cfg = []; cfg.outdir = save_dir; filename = ['OptTime_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

%% MEG vs. fMRI Correlation
close all
disp(nsub_IDs(discordantSubs)')
subjectForPlot = discordantSubs;
[optimalTimePoints, LI_opt] = findIndividualOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, subjectForPlot, lowerBound, upperBound);
cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
if LI_method == 1
    plotIndividualSourcePower(power_left, power_right, MEG_LI, fMRI_LI, sub_IDs, timePoints, LI_opt, subjectForPlot)
    cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_Dicordant_Pow_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
end

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

%% Response vs. Discordant LIs
discordant_subs = sub_MF_pt(discordantSubs);
discordant_indices = ismember(T_patn_MEGfMRI.Sub_ID, discordant_subs);
discordant_RT = T_patn_MEGfMRI.Avg(discordant_indices);
meanRT = nanmean(discordant_RT);

xllabel = [];
for i=1:length(discordantSubs)
    xllabel{i} = ['Subj ', num2str(discordantSubs(i))];
end

figure;
bar(discordant_RT);
xlabel('Subjects');
set(gca,'color','none');
ylabel('Response Time (sec)');
title({'Response Time', 'of Discordant Subjects'});
set(gca, 'XTick', 1:length(discordant_subs), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
set(gcf, 'Position', [1000   100   300   300]);
hold on;
line(get(gca,'xlim'), [meanRT meanRT], 'Color', 'red', 'LineStyle', '--');
hold off;
cfg = []; cfg.outdir = save_dir; filename = ['ReactionTime_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

%% Response Time Data (2)
figure;
daboxplot(T_patn_MEGfMRI.Avg, 'groups', ones(1, numel(T_patn_MEGfMRI.Avg)), 'outsymbol', 'kx', 'xtlabels', 'test', 'fill', 0);
xlabel('All Subjects');
ylabel('Response Time (sec)');
set(gca, 'FontSize', 10);
hold on;
hScatter = scatter(ones(size(discordant_RT)), discordant_RT, 'r', 'filled');
hold off;
title({'Response Time'});
l = legend(hScatter, 'Discordant', 'Location', 'south');
legendPos = l.Position;
legendPos(1) = legendPos(1) + 0.02;
l.Position = legendPos;
set(gca, 'XTick', []);
box off;
set(gcf, 'Position', [1000   100   300   300]);
cfg = []; cfg.outdir = save_dir; filename = ['2_ReactionTime_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

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
disp(['The total mean of mean_Falsefont_ACC is: ', num2str(totalmean)]);
totalstd = std(meanAccBySubject_Falsefont.mean_Falsefont_ACC);
disp(['The total std of mean_Falsefont_ACC is: ', num2str(totalstd)]);

%% Task Performance of Discordant Subjects
discordant_subs_numeric = cellfun(@(x) str2double(x(3:end)), discordant_subs);
discordant_indices_anim = find(ismember(meanAccBySubject_Animal.Subject, discordant_subs_numeric));
discordant_indices_symb = find(ismember(meanAccBySubject_Falsefont.Subject, discordant_subs_numeric));
meanTP_anim = nanmean(meanAccBySubject_Animal.mean_Animal_ACC);
meanTP_symb = nanmean(meanAccBySubject_Falsefont.mean_Falsefont_ACC);

% Plotting Animal Task Performance
figure;
bar(meanAccBySubject_Animal.mean_Animal_ACC(discordant_indices_anim));
xlabel('Subjects');
set(gca,'color','none');
ylabel('Accuracy (%)');
title({'Task performance (anim)', 'of Discordant Subjects'});
% set(gca, 'XTick', 1:length(discordant_indices_anim), 'XTickLabel', meanAccBySubject_Animal.Subject(discordant_indices_anim), 'XTickLabelRotation', 45);
set(gca, 'XTick', 1:length(discordant_indices_anim), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
set(gcf, 'Position', [1000   100   300   300]);
hold on;
line(get(gca,'xlim'), [meanTP_anim meanTP_anim], 'Color', 'red', 'LineStyle', '--');
hold off;
cfg = []; cfg.outdir = save_dir; filename = ['TaskPerformace_anim_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

% Plotting Symbol Task Performance
figure;
bar(meanAccBySubject_Falsefont.mean_Falsefont_ACC(discordant_indices_symb));
xlabel('Subjects');
set(gca,'color','none');
ylabel('Accuracy (%)');
title({'Task performance (fonts)', 'of Discordant Subjects'});
set(gca, 'XTick', 1:length(discordant_indices_symb), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
set(gcf, 'Position', [1000   100   300   300]);
hold on;
line(get(gca,'xlim'), [meanTP_symb meanTP_symb], 'Color', 'red', 'LineStyle', '--');
hold off;
cfg = []; cfg.outdir = save_dir; filename = ['TaskPerformace_symb_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

%% Task Performance (2)
figure;
daboxplot(meanAccBySubject_Animal.mean_Animal_ACC, 'groups', ones(1, numel(meanAccBySubject_Animal.mean_Animal_ACC)), 'outsymbol', 'kx', 'xtlabels', 'test', 'fill', 0);
xlabel('All Subjects');
ylabel('Accuracy (%)');
set(gca, 'FontSize', 10);
discordant_ACC_anim = meanAccBySubject_Animal.mean_Animal_ACC(discordant_indices_anim);
hold on;
hScatter = scatter(ones(size(discordant_ACC_anim)), discordant_ACC_anim, 'r', 'filled');
hold off;
title({'Task performance (anim)', 'of Discordant Subjects'});
l = legend(hScatter, 'Discordant', 'Location', 'south');
legendPos = l.Position;
legendPos(1) = legendPos(1) + 0.02;
l.Position = legendPos;
set(gca, 'XTick', []);
box off;
set(gcf, 'Position', [1000   100   300   300]);
set(gca,'color','none');
cfg = []; cfg.outdir = save_dir; filename = ['2_TaskPerformace_anim_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

figure;
daboxplot(meanAccBySubject_Falsefont.mean_Falsefont_ACC, 'groups', ones(1, numel(meanAccBySubject_Falsefont.mean_Falsefont_ACC)), 'outsymbol', 'kx', 'xtlabels', 'test', 'fill', 0);
xlabel('All Subjects');
ylabel('Accuracy (%)');
set(gca, 'FontSize', 10);
discordant_ACC_symb = meanAccBySubject_Falsefont.mean_Falsefont_ACC(discordant_indices_symb);
hold on;
hScatter = scatter(ones(size(discordant_ACC_symb)), discordant_ACC_symb, 'r', 'filled');
hold off;
title({'Task performance (symb)', 'of Discordant Subjects'});
l = legend(hScatter, 'Discordant', 'Location', 'south');
legendPos = l.Position;
legendPos(1) = legendPos(1) + 0.02;
l.Position = legendPos;
set(gca, 'XTick', []);
box off;
set(gcf, 'Position', [1000   100   300   300]);
set(gca,'color','none');

cfg = []; cfg.outdir = save_dir; filename = ['2_TaskPerformace_symb_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

%% Neuropsych Details of Discordant Patients
disp(patn_MEGfMRI_neuropsych(discordant_indices_anim,:))

%% Corr (Response ACC, Reaction Time)
% Convert numerical subject IDs in accuracyResults_updt to string format with 'EC' prefix
accuracyResults_updt.SubjectStr = strcat('EC', arrayfun(@num2str, accuracyResults_updt.Subject, 'UniformOutput', false));

% Find the intersection between subject IDs
[commonSubIDs, idxAccuracy, idxRT] = intersect(accuracyResults_updt.SubjectStr, T_patn_MEGfMRI.Sub_ID);

% Extract relevant data
animalAcc = accuracyResults_updt.Animal_ACC(idxAccuracy);
falsefontAcc = accuracyResults_updt.Falsefont_ACC(idxAccuracy);
animalRT = T_patn_MEGfMRI.Animal(idxRT);
falsefontRT = T_patn_MEGfMRI.Symbol(idxRT);

% Calculate correlation for animal condition
[correlationAnimal, pValueAnimal] = corr(animalAcc, animalRT, 'Rows', 'complete');

% Calculate correlation for falsefont condition
[correlationFalsefont, pValueFalsefont] = corr(falsefontAcc, falsefontRT, 'Rows', 'complete');

% Display results
disp(['Correlation between Animal Accuracy and Reaction Time: ', num2str(correlationAnimal)]);
disp(['P-value for Animal Accuracy and Reaction Time: ', num2str(pValueAnimal)]);

disp(['Correlation between Falsefont Accuracy and Reaction Time: ', num2str(correlationFalsefont)]);
disp(['P-value for Falsefont Accuracy and Reaction Time: ', num2str(pValueFalsefont)]);

% Plot Correlations

% Plot for Animal condition
figure;
subplot 211
scatter(animalAcc, animalRT, 'filled');
title('Corr (Animal Acc, ReactionT');
xlabel('Animal Acc');
ylabel('ReactionT (s)');
grid on;
set(gca,'color','none');

% Plot for Falsefont condition
subplot 212
scatter(falsefontAcc, falsefontRT, 'filled');
title('Corr (Falsefont Acc, ReactionT)');
xlabel('Falsefont Acc');
ylabel('ReactionT (s)');
grid on;
box off;
set(gcf, 'Position', [1000   100   300   300]);
set(gca,'color','none');

cfg = []; cfg.outdir = save_dir; filename = ['Corr_ReactionTvsAcc_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

cd(save_dir)

