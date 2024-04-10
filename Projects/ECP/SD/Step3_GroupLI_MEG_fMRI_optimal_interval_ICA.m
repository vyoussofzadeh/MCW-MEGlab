
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
    case 4
        cfg.datamask = fullfile('./Group_analysis/LCMV/results_abs*.mat');
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

%% avg.
% sFiles = S_data_anim_pt.sFiles_in;
% % Process: Average: Everything
% sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
%     'avgtype',         1, ...  % Everything
%     'avg_func',        1, ...  % Arithmetic average:  mean(x)
%     'weighted',        0, ...
%     'scalenormalized', 0);

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

%% MEG vs. fMRI lat analysis (PT)
% pause, close all,

[sub_MF_pt,IA,IB] = intersect(LI_pt_ID, fmri_LIs.ID.language_Lateral');

LI_pt_val = LI_pt.LI_sub;

MEG_LI_Data = LI_pt_val(:,IA,:);
fMRI_LI = (fmri_LIs.val.language_Lateral(IB)); % Lateral regions was used.

% missing from fMRI
difference = setdiff(LI_pt_ID, sub_MF_pt');
disp('missing from fMRI')
disp(difference');

%% MEG LI vs fMRI LI (language_Lateral)
% Corr, MEG-fMRI
cfg = []; cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.ternary = 0;
cfg.thre = .2;
cfg.savefig = 1;
cfg.bf = 5;
cfg.outdir = save_dir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = MEG_LI_Data;
cfg.fmri_LIs_val = fMRI_LI;
% cfg.net_sel = [1,2,6];
cfg.net_sel = [11]; % 6, 1, 2
[megLI_sub_pt, fmri_LIs_val, ~, interval_idx] = do_MEG_fMRI_corr_contrast(cfg);

%%
cfg = [];
cfg.thre = .2; cfg.LI = fmri_LIs_val;
fmri_LIs_trn = do_ternary_classification(cfg);
size(fmri_LIs_trn);

%%
% Define thresholds for MEG and fMRI
MEG_thre = 0.2; % MEG threshold
fMRI_thre = 0.1; % fMRI threshold

%%
clc
MEG_LI = squeeze(MEG_LI_Data(11,:,:)); MEG_LI(isnan(MEG_LI))=0;

%% fixed intervals LIs.
clc
timePoints = mean(wi,2);
IntervalSize = 1;
[optimalInterval, correlations] = findOptimalMEGInterval(MEG_LI, fMRI_LI, timePoints, IntervalSize);

% max(correlations)

sub_IDs = sub_MF_pt;
nsub_IDs = [];
for i=1:length(sub_IDs)
    nsub_IDs{i} = [num2str(i), ':', sub_IDs{i}];
end

optimalInterval_constant = ones (size(timePoints,1),1).*optimalInterval;

bf = 0.1;
lowerBound = optimalInterval - bf;
upperBound =  optimalInterval + bf;
[groupCorrelation_cnst, optimalTimePoints_cnst] = computeGroupLevelMEGfMRICorrelation_timepoints(MEG_LI, fMRI_LI, timePoints, lowerBound, upperBound);

[concordance_cnst, discordantSubs_cnst] = calculateConcordanceForTimePoints(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalInterval_constant);

% Display IDs of discordant subjects, if any
if ~isempty(discordantSubs_cnst)
    disp('Discordant Subjects:');
    disp(nsub_IDs(discordantSubs_cnst));
else
    disp('No Discordant Subjects Found');
end

disp(nsub_IDs');
nsub_IDs(discordantSubs_cnst)'
subjectForPlot = input('enter sub number:');
findIndividualOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, subjectForPlot, lowerBound, upperBound);
cfg = []; cfg.outdir = save_dir; filename = ['ConstantTimePoints_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

%%
network_sel = [1:3,6:11];
Data_hcp_atlas.groups_labels

%% Optimal interval LIs
close all
if LI_method == 1 % Power values
    
    Data_hcp_atlas.groups_labels

    pow_hc = transformPowSubTo3DArrays(LI_hc.pow_sub);
    pow_pt = transformPowSubTo3DArrays(LI_pt.pow_sub);
    
    mPow_sub_hc_left = squeeze(nanmean(pow_hc.left,2)); mPow_sub_hc_right = squeeze(nanmean(pow_hc.right,2));
    mPow_sub_pt_left = squeeze(nanmean(pow_pt.left,2)); mPow_sub_pt_right = squeeze(nanmean(pow_pt.right,2));
    
    power_left = squeeze(pow_pt.left(11,IA,:));
    power_right = squeeze(pow_pt.right(11,IA,:));

    fMRI_LI = (fmri_LIs.val.language_Lateral(IB)); % Lateral regions was used.

    %     close all
    plot_flag = 1;
    % Define the time interval bounds
    lowerBound = 0.5; %
    upperBound = 1.0; %
    %     optimalTimePoints_pow = plotSubjectPowerOverlay(power_left, power_right, sub_IDs, timePoints, plot_flag);
    [optimalTimePoints_pow] = plotSubjectPowerOverlay(power_left, power_right, sub_IDs, timePoints, plot_flag, lowerBound, upperBound);
    
    if plot_flag == 1
        cfg = []; cfg.outdir = save_dir; filename = ['Power_values_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig';
        do_export_fig(cfg)
    end
    
    [correlationCoeff, pValue, optimalMEGLI, optimalfMRILI] = calculateCorrelationForOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, optimalTimePoints_pow);
    [groupCorrelation, optimalTimePoints] = computeGroupLevelMEGfMRICorrelation_timepoints(MEG_LI, fMRI_LI, timePoints, lowerBound, upperBound);
    % Calculate concordance for the identified time points and identify discordant subjects
    [concordance_pow, discordantSubs_pow] = calculateConcordanceForTimePoints(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimePoints_pow);
    
    % Display concordance rate
    disp(['Concordance: ', num2str(discordantSubs_pow)]);
    disp(['Corr: ', num2str(groupCorrelation)]);

    
    % Display IDs of discordant subjects, if any
    if ~isempty(discordantSubs_pow)
        disp('Discordant Subjects:');
        disp(nsub_IDs(discordantSubs_pow));
    else
        disp('No Discordant Subjects Found');
    end
    
    disp(nsub_IDs');
    nsub_IDs(discordantSubs_pow)'
    subjectForPlot = input('enter sub number:');
    findIndividualOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, subjectForPlot, lowerBound, upperBound);
    cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    plotIndividualSourcePower(power_left, power_right, MEG_LI, sub_IDs, timePoints, subjectForPlot)
    cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_Dicordant_LIs_pow_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    
    % Number of subjects
    numSubjects = length(optimalTimePoints_pow);
    
    % Generate subject IDs (if not already provided)
    subjectIDs = 1:numSubjects; % Adjust as necessary based on your data structure
    
    % Create a scatter plot of optimal time points for each subject
    figure; scatter(subjectIDs, optimalTimePoints_pow, 'filled');
    title('Optimal Time Points for Each Subject');
    xlabel('Subject ID');
    ylabel('Optimal Time Point (s)');
    grid on; % Add a grid for easier visualization
    set(gca,'color','none');
    
    % L = length(sub_MF_pt);
    % set(gca,'Xtick', 1:L,'XtickLabel',sub_MF_pt);
    % set(gca,'FontSize',8,'XTickLabelRotation',90);
    
    cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    
    % Response time data
    % source code: /data/MEG/Research/aizadi/Scripts/Pipelines/ReactionTimeAnalysis/Pipe_process_task_RT_summary.m
    load('/data/MEG/Research/aizadi/process/RT_summary/ResponseTime.mat')
    
    discordant_subs = sub_MF_pt(discordantSubs_pow);
    
    % Find indices of discordant subjects in RT data
    discordant_indices = find(ismember(T.Sub_ID, discordant_subs));
    
    % Extract RT data for discordant subjects
    discordant_RT = T.Avg(discordant_indices);
    
    % Calculate the mean response time
    meanRT = nanmean(discordant_RT);
    
    % Plotting
    figure;
    bar(discordant_RT);
    xlabel('Subjects');
    set(gca,'color','none');
    ylabel('Response Time (sec)');
    title('Response Time of Discordant Subjects');
    set(gca, 'XTick', 1:length(discordant_subs), 'XTickLabel', discordant_subs, 'XTickLabelRotation', 45);
    set(gcf, 'Position', [1000   100   300   300]);
    
    % Draw a horizontal line at the mean response time
    hold on; % Keep the current bar plot
    line(get(gca,'xlim'), [meanRT meanRT], 'Color', 'red', 'LineStyle', '--'); % Draw line
    hold off; % Release the plot
    
    cfg = []; cfg.outdir = save_dir; filename = ['ReactionTime_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    cd(save_dir)
    
elseif LI_method == 2 % counting

    Data_hcp_atlas.groups_labels
    fMRI_LI = (fmri_LIs.val.language_Lateral(IB)); % Lateral regions was used.
   
    roicnt_hc = (LI_hc.roi_count_sub);
    roicnt_pt = (LI_pt.roi_count_sub);
    
    roicnt_pt_lat = squeeze(roicnt_pt(11,IA,:));
    
    %     close all
    plot_flag = 0;
    % Define the time interval bounds
    lowerBound = 0.5; %
    upperBound = 1.0; %
    [optimalTimePoints] = plotSubjectActiveverciesOverlay(roicnt_pt_lat, sub_IDs, timePoints, plot_flag, lowerBound, upperBound);
    
    
    if plot_flag == 1
        cfg = []; cfg.outdir = save_dir; filename = ['Counting_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig';
        do_export_fig(cfg)
    end
    
    [correlationCoeff, pValue, optimalMEGLI, optimalfMRILI] = calculateCorrelationForOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, optimalTimePoints);
    [groupCorrelation, optimalTimePoints] = computeGroupLevelMEGfMRICorrelation_timepoints(MEG_LI, fMRI_LI, timePoints, lowerBound, upperBound);


    % Calculate concordance for the identified time points and identify discordant subjects
    [concordance_count, discordantSubs_count] = calculateConcordanceForTimePoints(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimePoints);
    
    % Display concordance rate
    disp(['Concordance: ', num2str(discordantSubs_count)]);
    disp(['Corr: ', num2str(groupCorrelation)]);

    
    % Display IDs of discordant subjects, if any
    if ~isempty(discordantSubs_count)
        disp('Discordant Subjects:');
        disp(nsub_IDs(discordantSubs_count));
    else
        disp('No Discordant Subjects Found');
    end
    
    disp(nsub_IDs');
    nsub_IDs(discordantSubs_count)'
    subjectForPlot = input('enter sub number:');
    findIndividualOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, subjectForPlot, lowerBound, upperBound);
    cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    plotIndividualCounting(roicnt_pt_lat, MEG_LI, sub_IDs, timePoints, subjectForPlot)
    cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_Dicordant_LIs_count_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    
    % Number of subjects
    numSubjects = length(optimalTimePoints);
    
    % Generate subject IDs (if not already provided)
    subjectIDs = 1:numSubjects; % Adjust as necessary based on your data structure
    
    % Create a scatter plot of optimal time points for each subject
    figure; scatter(subjectIDs, optimalTimePoints, 'filled');
    title('Optimal Time Points for Each Subject');
    xlabel('Subject ID');
    ylabel('Optimal Time Point (s)');
    grid on; % Add a grid for easier visualization
    set(gca,'color','none');
    
    % L = length(sub_MF_pt);
    % set(gca,'Xtick', 1:L,'XtickLabel',sub_MF_pt);
    % set(gca,'FontSize',8,'XTickLabelRotation',90);
    
    cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)  
    
    % Response time data
    % source code: /data/MEG/Research/aizadi/Scripts/Pipelines/ReactionTimeAnalysis/Pipe_process_task_RT_summary.m
    load('/data/MEG/Research/aizadi/process/RT_summary/ResponseTime.mat')
    
    discordant_subs = sub_MF_pt(discordantSubs_count);
    
    % Find indices of discordant subjects in RT data
    discordant_indices = find(ismember(T.Sub_ID, discordant_subs));
    
    % Extract RT data for discordant subjects
    discordant_RT = T.Avg(discordant_indices);
    
    % Calculate the mean response time
    meanRT = nanmean(discordant_RT);
    
    % Plotting
    figure;
    bar(discordant_RT);
    xlabel('Subjects');
    set(gca,'color','none');
    ylabel('Response Time (sec)');
    title('Response Time of Discordant Subjects');
    set(gca, 'XTick', 1:length(discordant_subs), 'XTickLabel', discordant_subs, 'XTickLabelRotation', 45);
    set(gcf, 'Position', [1000   100   300   300]);
    
    % Draw a horizontal line at the mean response time
    hold on; % Keep the current bar plot
    line(get(gca,'xlim'), [meanRT meanRT], 'Color', 'red', 'LineStyle', '--'); % Draw line
    hold off; % Release the plot
    
    cfg = []; cfg.outdir = save_dir; filename = ['ReactionTime_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    cd(save_dir)
    
elseif LI_method == 3 % bootstrap

    Data_hcp_atlas.groups_labels
    fMRI_LI = (fmri_LIs.val.language_Temporal(IB)); % Lateral regions was used.
    
    MEG_LI = squeeze(MEG_LI_Data(6,:,:)); 
    
    
    % Optimal time points LIs.
    % Clear the command window and close all figures
    disp('========')
    clc
    % close all
    
    % Define the time interval bounds
    lowerBound = 0.5; %
    upperBound = 1; %
    
    % Compute group-level MEG-fMRI correlation and find optimal time points
    [groupCorrelation, optimalTimePoints] = computeGroupLevelMEGfMRICorrelation_timepoints(MEG_LI, fMRI_LI, timePoints, lowerBound, upperBound);
    
    % Plot optimal time points on MEG data
    plotOptimalTimePointsOnMEG(MEG_LI, fMRI_LI, sub_IDs, timePoints, optimalTimePoints);
    
    % Calculate and display the mean of the optimal time points
    meanOptimalTime = mean(optimalTimePoints);
    disp(['Mean of optimal time points: ', num2str(meanOptimalTime)]);
    disp(['Corr: ', num2str(groupCorrelation)]);

%     [correlationCoeff, pValue, optimalMEGLI, optimalfMRILI] = calculateCorrelationForOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, optimalTimePoints);
    
    % Calculate concordance for the identified time points and identify discordant subjects
    [concordance, discordantSubs] = calculateConcordanceForTimePoints(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimePoints);
    
    % Display concordance rate
    disp(['Concordance: ', num2str(concordance)]);
    
    % Display IDs of discordant subjects, if any
    if ~isempty(discordantSubs)
        disp('Discordant Subjects:');
        disp(nsub_IDs(discordantSubs));
    else
        disp('No Discordant Subjects Found');
    end
    
    disp(nsub_IDs');
    nsub_IDs(discordantSubs)'
    subjectForPlot = input('enter sub number:');
    findIndividualOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, subjectForPlot, lowerBound, upperBound);
    cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    % Number of subjects
    numSubjects = length(optimalTimePoints);
    
    % Generate subject IDs (if not already provided)
    subjectIDs = 1:numSubjects; % Adjust as necessary based on your data structure
    
    % Create a scatter plot of optimal time points for each subject
    figure; scatter(subjectIDs, optimalTimePoints, 'filled');
    title('Optimal Time Points for Each Subject');
    xlabel('Subject ID');
    ylabel('Optimal Time Point (s)');
    grid on; % Add a grid for easier visualization
    set(gca,'color','none');
    
    % L = length(sub_MF_pt);
    % set(gca,'Xtick', 1:L,'XtickLabel',sub_MF_pt);
    % set(gca,'FontSize',8,'XTickLabelRotation',90);
    
    cfg = []; cfg.outdir = save_dir; filename = ['optimalTimePoints_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    
    % Response time data
    % source code: /data/MEG/Research/aizadi/Scripts/Pipelines/ReactionTimeAnalysis/Pipe_process_task_RT_summary.m
    load('/data/MEG/Research/aizadi/process/RT_summary/ResponseTime.mat')
    
    discordant_subs = sub_MF_pt(discordantSubs);
    
    % Find indices of discordant subjects in RT data
    discordant_indices = find(ismember(T.Sub_ID, discordant_subs));
    
    % Extract RT data for discordant subjects
    discordant_RT = T.Avg(discordant_indices);
    
    % Calculate the mean response time
    meanRT = nanmean(discordant_RT);
    
    % Plotting
    figure;
    bar(discordant_RT);
    xlabel('Subjects');
    set(gca,'color','none');
    ylabel('Response Time (sec)');
    title('Response Time of Discordant Subjects');
    set(gca, 'XTick', 1:length(discordant_subs), 'XTickLabel', discordant_subs, 'XTickLabelRotation', 45);
    set(gcf, 'Position', [1000   100   300   300]);
    
    % Draw a horizontal line at the mean response time
    hold on; % Keep the current bar plot
    line(get(gca,'xlim'), [meanRT meanRT], 'Color', 'red', 'LineStyle', '--'); % Draw line
    hold off; % Release the plot
    
    cfg = []; cfg.outdir = save_dir; filename = ['ReactionTime_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    cd(save_dir)
end
