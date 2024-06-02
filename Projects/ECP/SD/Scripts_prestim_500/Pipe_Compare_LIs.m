
%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 12/06/2023

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/run')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/function')
Run_setpath

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other/')

%%
% LI_analysis_label = {'DICS_baseline','DICS_contrast','LCMV_baseline','LCMV_contrast','DICS_anim', 'DICS_contrast_prestim', 'dSPM_contrast'};
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
% data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results';
data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim';
cd(data_save_dir)

%%
% save_dir = fullfile(data_save_dir,LI_analysis_label{LI_analysis}, 'compare_LIs');
% checkOrCreateDir(save_dir)
% cd(save_dir)

save_dir = fullfile(data_save_dir,LI_analysis_label{LI_analysis}, 'compare_LIs');
checkOrCreateDir(save_dir)
cd(save_dir)

LI_pt = [];
for i=1:length(LI_method_label)
    LI_pt.(LI_method_label{i}) = load(fullfile(data_save_dir,LI_analysis_label{LI_analysis},LI_method_label{i},'LI_Patn'));
end

%% Run Brainstorm
% Run_BS
% %% BS
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
        %     case 4
        %         cfg.datamask = fullfile('./Group_analysis/LCMV/results_abs*.mat');
        %         S_data = ecpfunc_read_sourcemaps_contrast(cfg);
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
    case {2,4,6}
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
[sub_MF_pt,IA,IB] = intersect(LI_pt_ID, fmri_LIs.ID.language_Lateral');

LI_pt_val = LI_pt.(LI_method_label{1}).LI_sub;

% LI_pt_val_new = LI_pt_val(:,IA,:);
fmri_LIs_val = (fmri_LIs.val.language_Lateral(IB)); % Lateral regions was used.

% missing from fMRI
difference = setdiff(LI_pt_ID, sub_MF_pt');
disp('missing from fMRI')
disp(difference');

%% MEG LI vs fMRI LI (language_Lateral)
LI_pt_val_new = [];
for i=1:length(LI_method_label)
    LI_pt_val_new.(LI_method_label{i}) = LI_pt.(LI_method_label{i}).LI_sub(:, IA,:);
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
[megLI_sub_pt, fmri_LIs_val, ~] = do_MEG_fMRI_corr_contrast_approches(cfg);
cfg.net_sel = [1];
[megLI_sub_pt, fmri_LIs_val, ~] = do_MEG_fMRI_corr_contrast_approches(cfg);
cfg.net_sel = [5];
[megLI_sub_pt, fmri_LIs_val, ~] = do_MEG_fMRI_corr_contrast_approches(cfg);
cfg.net_sel = [2];
[megLI_sub_pt, fmri_LIs_val, ~] = do_MEG_fMRI_corr_contrast_approches(cfg);
cfg.net_sel = [6];
[megLI_sub_pt, fmri_LIs_val, ~] = do_MEG_fMRI_corr_contrast_approches(cfg);
cfg.net_sel = [2,6];
[megLI_sub_pt, fmri_LIs_val, ~] = do_MEG_fMRI_corr_contrast_approches(cfg);
cfg.net_sel = [1,2,6];
[megLI_sub_pt, fmri_LIs_val, ~] = do_MEG_fMRI_corr_contrast_approches(cfg);
cfg.net_sel = [9];
[megLI_sub_pt, fmri_LIs_val, ~] = do_MEG_fMRI_corr_contrast_approches(cfg);

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
cfg.buffervalue = 1;
[megLIs_trn, fmri_LIs_trn] = do_MEG_fMRI_concordance_contrast_approches(cfg);
cfg.net_sel = [1];
[megLIs_trn, fmri_LIs_trn] = do_MEG_fMRI_concordance_contrast_approches(cfg);
cfg.net_sel = [5];
[megLIs_trn, fmri_LIs_trn] = do_MEG_fMRI_concordance_contrast_approches(cfg);
cfg.net_sel = [2];
[megLIs_trn, fmri_LIs_trn] = do_MEG_fMRI_concordance_contrast_approches(cfg);
cfg.net_sel = [6];
[megLIs_trn, fmri_LIs_trn] = do_MEG_fMRI_concordance_contrast_approches(cfg);
cfg.net_sel = [2,6];
[megLIs_trn, fmri_LIs_trn] = do_MEG_fMRI_concordance_contrast_approches(cfg);
cfg.net_sel = [1,2,6];
[megLIs_trn, fmri_LIs_trn] = do_MEG_fMRI_concordance_contrast_approches(cfg);
cfg.net_sel = [9];
[megLIs_trn, fmri_LIs_trn] = do_MEG_fMRI_concordance_contrast_approches(cfg);

%% mean MEG li vs. fMRI
% pause,
%
close all,
clc

for i=1:length(LI_method_label)
    
    disp(['correlation analysis: ',LI_method_label{i}])
    
    cfg = []; cfg.wi = wi;
    cfg.ID = sub_MF_pt;
    cfg.thre = 0.10;
    cfg.bf = 1;
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
    do_MEG_fMRI_corr_contrast_rois(cfg);
    
    % - export figs
    cfg = []; cfg.outdir = save_dir; filename = ['corr ROIs_', LI_method_label{i}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
    disp(['Concordance analysis: ',LI_method_label{i}])
    
    % clc
    cfg = []; cfg.wi = wi;
    cfg.ID = sub_MF_pt;
    cfg.thre = 0.10;
    cfg.ternary = 1;
    cfg.savefig = 0;
    cfg.outdir = save_dir;
    cfg.net_sel_mutiple_label = net_sel_mutiple_label;
    cfg.net_sel_id = [1,2,5,6,11];
    cfg.lang_id = {'language_Angular'; 'language_Frontal';'language_PCingPrecun'; 'language_Temporal'; 'language_Lateral'};
    cfg.LI_val = LI_pt_val_new.(LI_method_label{i});
    cfg.fmri_LIs_val = fmri_LIs_trn;
    cfg.idx = IB;
    cfg.title = LI_method_label{i};
    cfg.buffervalue = 5;
    do_MEG_fMRI_concordance_contrast_rois(cfg);
    
    % - export figs
    cfg = []; cfg.outdir = save_dir; filename = ['concor ROIs_', LI_method_label{i}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
    
end
cd(save_dir)

%%
close all
% Initialize a table to store results
resultsTable = table([], [], 'VariableNames', {'Method', 'Metrics'});

for i=1:length(LI_method_label)
    
    disp(['Processing: ', LI_method_label{i}]);
    
    % Configuration for both analyses
    cfg = [];
    cfg.wi = wi;
    cfg.ID = sub_MF_pt;
    cfg.savefig = 1; % Assuming this controls figure saving inside the functions
    cfg.outdir = save_dir;
    cfg.net_sel_mutiple_label = net_sel_mutiple_label;
    cfg.net_sel_id = [1,2,5,6,11];
    cfg.lang_id = {'language_Angular'; 'language_Frontal';'language_PCingPrecun'; 'language_Temporal'; 'language_Lateral'};
    cfg.LI_val = LI_pt_val_new.(LI_method_label{i});
    cfg.fmri_LIs_val = fmri_LIs;
    cfg.idx = IB;
    cfg.title = LI_method_label{i};
    
    % Correlation Analysis
    cfg.thre = 0.10;
    cfg.ternary = 0;
    [~, ~, correlationMetrics] = do_MEG_fMRI_corr_contrast_rois(cfg);
    
    % Concordance Analysis
    cfg.thre = 0.15;
    cfg.ternary = 1;
    cfg.buffervalue = 2;
    cfg.fmri_LIs_val = fmri_LIs_trn;
    [~, ~, concordanceMetrics] = do_MEG_fMRI_concordance_contrast_rois(cfg);
    
    % Compile results
    metrics = struct('Correlation', correlationMetrics, 'Concordance', concordanceMetrics);
    resultsTable = [resultsTable; {LI_method_label{i}, metrics}];
end

%%
close all

% Example metric name - adjust according to your actual data structure
metricName = {'Correlation','Concordance'}; % Assuming this is the name of the correlation/concordance field
roi_label = {'Angular'; 'Frontal';'PCingPrecun'; 'Temporal'; 'Lateral'};

% Colors for each method, adjust or extend as needed
colors = lines(length(LI_method_label));

% Loop through each method and plot
for k = 1:length(metricName)
    for j = 1:length(roi_label)
        % Initialize a figure
        figure;
        hold on; % Hold on to plot multiple lines
        for i = 1:height(resultsTable)
            % Extract the relevant metric for the current method
            currentMetrics = resultsTable.Metrics(i).(metricName{k});            
            % Example plotting command - adjust as needed
            plot(mean(wi'), currentMetrics(j,:), 'LineWidth', 3, 'Color', colors(i,:));
            title([resultsTable.Method{i}, '-', roi_label{j}]);
        end        
        
        % Beautify the plot
        set(gca, 'color', 'none');
        ylabel(['LIs ', metricName{k}]);
        xlabel('Time (sec)');
        legend(resultsTable.Method, 'Location', 'southoutside', 'NumColumns', 5);
        box off;
        hold off; % Release the hold to stop plotting on the same figure
        set(gcf, 'Position', [100   400   300  400]);
    end
end
% Optionally, set the title based on `net_sel_mutiple_label` and `net_sel_id`
% title([net_sel_mutiple_label{net_sel}]);

%%
clc
close all;

% Assuming metricName contains the names of the fields we want to plot
metricNames = {'Correlation', 'Concordance'};
roi_labels = {'Angular', 'Frontal', 'PCingPrecun', 'Temporal', 'Lateral'};

baseColors = [
    31, 78, 121;
    132, 60, 12;
    127, 96, 0;
    112, 48, 160;
    84, 130, 53]/ 255;

% Enhance color contrast
contrastFactorMagnitude = 0.5;  % Higher factor for lighter colors in Magnitude
contrastFactorCounting = 0.2;   % Lower factor for slightly brighter colors in Counting

colorsMagnitude = min(baseColors + contrastFactorMagnitude, 1);
colorsCounting = min(baseColors + contrastFactorCounting, 1);


% Loop through each metric type (Correlation and Concordance)
for metricIdx = 1:length(metricNames)
    figure;
    
    for roiIdx = 1:length(roi_labels)
        subplot(length(roi_labels), 1, roiIdx);
        hold on;
        
        maxValue = -inf;
        maxTime = 0;
        
        for i = 1:height(resultsTable)
            if isfield(resultsTable.Metrics(i), metricNames{metricIdx})
                currentMetrics = resultsTable.Metrics(i).(metricNames{metricIdx});
                
                if size(currentMetrics, 1) >= roiIdx
                    dataToPlot = currentMetrics(roiIdx, :);
                    
                    % Using if-else to determine the color to usedi
                    if i == 1
                        colorToUse = colorsMagnitude(roiIdx, :);
                    elseif i == 2
                        colorToUse = colorsCounting(roiIdx, :);
                    else
                        colorToUse = baseColors(roiIdx, :);
                    end
                    
                    plot(mean(wi'), dataToPlot, 'LineWidth', 2, 'Color', colorToUse);
                    
                    [localMax, localMaxIndex] = max(dataToPlot);
                    if localMax > maxValue
                        maxValue = localMax;
                        maxTime = mean(wi(localMaxIndex));
                    end
                end
            end
        end       
        
        if maxValue > -inf
            text(maxTime, maxValue, sprintf('Max: %.2f', maxValue), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        end
        
        set(gca, 'color', 'none');
        if strcmp(metricNames{metricIdx}, 'Correlation')
            ylim([-0.2 0.8]);
        elseif strcmp(metricNames{metricIdx}, 'Concordance')
            ylim([30 90]);
        end
        
        if roiIdx == length(roi_labels)
            ylabel(['LIs ', metricNames{metricIdx}, ' (MEG vs. fMRI)']);
            xlabel('Time (sec)');
        end
        legend(resultsTable.Method, 'Location', 'southoutside', 'NumColumns', 2, 'Orientation', 'horizontal');

        
        title(roi_labels{roiIdx});
        box off;
%         axis square
        hold off;
    end
    
    % Super title for the entire figure
    sgtitle([metricNames{metricIdx}, ' Analysis']);
    set(gcf, 'Position', [100, 400, 350, 800]);
    
    cfg = []; cfg.outdir = save_dir; filename = [metricNames{metricIdx}, ' rois']; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg);
end
disp(['saved as, ', filename])

%%
close all;

% Assuming metricName contains the names of the fields we want to plot
metricNames = {'Correlation', 'Concordance'};
roi_labels = {'Angular', 'Frontal', 'PCingPrecun', 'Temporal', 'Lateral'};

% Colors for each method, adjust or extend as needed
colors = lines(height(resultsTable));

% Loop through each ROI to create a figure
for roiIdx = 1:length(roi_labels)
    figure; % Initialize a figure for the current ROI
    sgtitle([roi_labels{roiIdx}]); % Super title for the entire figure
    set(gcf, 'Position', [100, 100, 400, 600]); % Adjust figure size
    
    % Create a subplot for each metric type within the current ROI figure
    for metricIdx = 1:length(metricNames)
        subplot(1, length(metricNames), metricIdx);
        hold on; % Hold on to plot multiple lines in the subplot
        
        % Plotting loop for each method
        for i = 1:height(resultsTable)
            % Check if the current metric exists in the current method's metrics
            if isfield(resultsTable.Metrics(i), metricNames{metricIdx})
                currentMetrics = resultsTable.Metrics(i).(metricNames{metricIdx});
                
                % Check if we have enough data for the current ROI
                if size(currentMetrics, 1) >= roiIdx
                    plot(mean(wi'), currentMetrics(roiIdx, :), 'LineWidth', 3, 'Color', colors(i,:));
                end
            end
        end
        
        % Plot adjustments for the current subplot
        set(gca, 'color', 'none');
        if strcmp(metricNames{metricIdx}, 'Correlation')
            ylim([-.2, 1]); % Set ylim for Correlation analysis
        elseif strcmp(metricNames{metricIdx}, 'Concordance')
            ylim([0, 90]); % Set ylim for Concordance analysis
        end
        
        % Apply labels only on relevant subplots
        if metricIdx == 1
            ylabel(['LIs ', '(MEG vs. fMRI)']);
        end
        xlabel('Time (sec)');
        
        % Add legend below the last subplot
        if roiIdx == length(roi_labels) && metricIdx == length(metricNames)
            lgd = legend(resultsTable.Method, 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 2);
            lgdPos = lgd.Position; % Get current position
            lgdPos(2) = lgdPos(2) - 0.11; % Move legend down
            lgd.Position = lgdPos; % Set new position
        end
        
        title(metricNames{metricIdx});
        box off;
        hold off; % Release the hold for the next metric
    end
end

%%
close all;

% Assuming metricName contains the names of the fields we want to plot
metricNames = {'Correlation', 'Concordance'};
roi_labels = {'Ang', 'Front', 'PCingPrecun', 'Temp', 'Lat'};
LI_method_labels = {'SourceMag', 'Count', 'Bootstrp'}; % LI Methods

% Colors for each ROI, adjust or extend as needed
colors = lines(length(roi_labels));

% Loop through each LI method
for methodIdx = 1:length(LI_method_labels)
    figure; % Initialize a figure for the current LI method
    sgtitle([LI_method_labels{methodIdx}, ' analysis: ROIs']); % Super title for the figure
    set(gcf, 'Position', [100, 100, 600, 600]); % Adjust figure size
    
    % Create subplots for Correlation and Concordance
    for metricIdx = 1:length(metricNames)
        subplot(1, 2, metricIdx);
        hold on; % Hold on to overlay multiple lines in the subplot
        
        % Loop through each ROI to overlay in the current subplot
        for roiIdx = 1:length(roi_labels)
            % Assuming there's a way to access the correct set of metrics for the current LI method
            currentMetrics = []; % Initialize empty; this needs to be filled with actual data retrieval logic
            if isfield(resultsTable.Metrics(methodIdx), metricNames{metricIdx})
                currentMetrics = resultsTable.Metrics(methodIdx).(metricNames{metricIdx});
                
                % Check if we have enough data for the current ROI
                if size(currentMetrics, 1) >= roiIdx
                    plot(mean(wi'), currentMetrics(roiIdx, :), 'LineWidth', 3, 'Color', colors(roiIdx,:));
                end
            end
        end
        
        % Plot adjustments for the current subplot
        set(gca, 'color', 'none');
        xlabel('Time (sec)');
        
        % Apply specific y-axis limits based on the metric
        if strcmp(metricNames{metricIdx}, 'Correlation')
            ylim([-0.3, 1]); % Set ylim for Correlation analysis
        elseif strcmp(metricNames{metricIdx}, 'Concordance')
            ylim([0, 90]); % Set ylim for Concordance analysis
        end
        
        title(metricNames{metricIdx});
        box off;
        
        if metricIdx == length(metricNames) % Add legend only to the last subplot
            lgd = legend(roi_labels, 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', length(roi_labels));
            lgdPos = lgd.Position; % Get current position
            lgdPos(2) = lgdPos(2) - 0.10; % Move legend down
            lgd.Position = lgdPos; % Set new position
        end
        
        hold off; % Release hold for next metric type subplot
    end
    cfg = []; cfg.outdir = save_dir; filename = [LI_method_labels{methodIdx}, ' Analysis across ROIs']; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)
end

%%
% Initialize a table to store the summary
summaryTable = table();

% Loop through LI methods, metric types, and ROIs
for methodIdx = 1:length(LI_method_labels)
    for metricIdx = 1:length(metricNames)
        for roiIdx = 1:length(roi_labels)
            currentMetrics = resultsTable.Metrics(methodIdx).(metricNames{metricIdx})(roiIdx, :);
            [maxValue, maxIndex] = max(currentMetrics); % Find max value and its index
            
            % Calculate corresponding time from 'wi' using 'maxIndex'
            maxTime = mean(wi(maxIndex, :), 2); % Assuming 'wi' contains start and end times of intervals
            
            % Add to the summary table
            newRow = {LI_method_labels{methodIdx}, metricNames{metricIdx}, roi_labels{roiIdx}, maxValue, maxTime};
            summaryTable = [summaryTable; newRow];
        end
    end
end

% Set column names
summaryTable.Properties.VariableNames = {'LI_Method', 'Metric_Type', 'ROI', 'Max_Value', 'Time_Interval'};

writetable(summaryTable, 'LI_Metrics_Summary.csv');

%%
fid = fopen('LI_Metrics_Summary.txt', 'wt'); % Open file for writing
% Print a header
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', summaryTable.Properties.VariableNames{:});

% Loop through each row of the table and print
for i = 1:height(summaryTable)
    fprintf(fid, '%s\t%s\t%s\t%f\t%f\n', summaryTable.LI_Method{i}, summaryTable.Metric_Type{i}, ...
        summaryTable.ROI{i}, summaryTable.Max_Value(i), summaryTable.Time_Interval(i));
end

fclose(fid); % Close the file


%%
close all
% Unique ROIs for iteration
uniqueROIs = unique(summaryTable.ROI);

% Loop through each ROI for plotting
figure;
sgtitle('Corr Max'); % Super title for the figure
for i = 1:length(uniqueROIs)
    
    subplot (3,2,i)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = summaryTable(strcmp(summaryTable.ROI, roi), :);
    
    % Create figure for current ROI
    
    k=1;
    % Plot Max Values for Correlation
    hold on;
    for method = unique(roiData.LI_Method)'
        methodData = roiData(strcmp(roiData.LI_Method, method) & strcmp(roiData.Metric_Type, 'Correlation'), :);
        bar(categorical(methodData.LI_Method), methodData.Max_Value, 'BarWidth', 0.2);
        text(k, methodData.Max_Value + 0.02, sprintf('%.1f', methodData.Max_Value), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
        k=1+k;
    end
    title([roi]);
    ylabel('Corr.');
    hold off;
    set(gca,'color','none');
    axis tight
    ylim([0, 1])
    set(gcf, 'Position', [100, 100, 400, 800]); % Adjust figure size
end

cfg = []; cfg.outdir = save_dir; filename = ' Corr max'; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)


figure;
sgtitle('Conc Max'); % Super title for the figure
for i = 1:length(uniqueROIs)
    
    subplot (3,2,i)
    % Plot Max Values for Concordance
    
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = summaryTable(strcmp(summaryTable.ROI, roi), :);
    
    % Create figure for current ROI
    k=1;
    hold on;
    for method = unique(roiData.LI_Method)'
        methodData = roiData(strcmp(roiData.LI_Method, method) & strcmp(roiData.Metric_Type, 'Concordance'), :);
        bar(categorical(methodData.LI_Method), methodData.Max_Value, 'BarWidth', 0.2);
        text(k, methodData.Max_Value + 0.02, sprintf('%.1f', methodData.Max_Value), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
        k=1+k;
    end
    title([roi]);
    ylabel('Concordance');
    hold off;
    set(gca,'color','none');
    ylim([0, 100])
    set(gcf, 'Position', [100, 100, 400, 800]); % Adjust figure size
    
end

cfg = []; cfg.outdir = save_dir; filename = ' Conc max'; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)


%%
close all
% Loop through each ROI for plotting
figure;
sgtitle('Correlation Max Values'); % Super title for the figure
for i = 1:length(uniqueROIs)
    subplot(3, 2, i);
    roi = uniqueROIs{i};

    % Extract data for the current ROI for Correlation
    roiData = summaryTable(strcmp(summaryTable.ROI, roi) & strcmp(summaryTable.Metric_Type, 'Correlation'), :);
    
    % Assuming roiData.LI_Method is categorical with unique categories for each method
    % If not, consider converting or ensuring uniqueness
    methods = categories(categorical(roiData.LI_Method));
    maxValues = roiData.Max_Value;
    
    b = bar(1:length(methods), maxValues, 'BarWidth', 0.2); % Plot bars
    ytickangle(45); xtickangle(45);
    
    % Annotations
    for j = 1:length(maxValues)
        % Manually calculate X position (center of each bar) and Y position (height of each bar)
        x = b.XData(j);
        y = maxValues(j);
        text(x, y, sprintf('%.2f', y), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
    
    set(gca, 'xtick', 1:length(methods), 'xticklabel', methods);
    title(roi);
    ylabel('Corr.');
    set(gca, 'color', 'none');
    ylim([0, 1]); % Adjust based on your data range
end
set(gcf, 'Position', [100, 100, 400, 700]); % Adjust figure size

cfg = []; cfg.outdir = save_dir; filename = ' Corr max rois'; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

figure;
sgtitle('Concordance Max Values'); % Super title for the figure
for i = 1:length(uniqueROIs)
    subplot(3, 2, i);
    roi = uniqueROIs{i};

    % Extract data for the current ROI for Concordance
    roiData = summaryTable(strcmp(summaryTable.ROI, roi) & strcmp(summaryTable.Metric_Type, 'Concordance'), :);
    
    % Ensure LI_Method is categorical and methods are unique for plotting
    methods = categories(categorical(roiData.LI_Method));
    maxValues = roiData.Max_Value; % Extract max values for Concordance
    
    b = bar(1:length(methods), maxValues, 'BarWidth', 0.4); % Plot bars for Concordance
    ytickangle(45); xtickangle(45);
    
    % Annotations
    for j = 1:length(maxValues)
        x = b.XData(j); % Use XData for horizontal position
        y = maxValues(j); % Max value for vertical position
       text(x, y + 2, sprintf('%.1f', y), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
    
    set(gca, 'xtick', 1:length(methods), 'xticklabel', methods);
    title([roi]);
    ylabel('Concordance');
    set(gca, 'color', 'none');
    ylim([0, 100]); % Adjust based on your range of Concordance values
end
set(gcf, 'Position', [100, 100, 400, 700]); % Adjust figure size

cfg = []; cfg.outdir = save_dir; filename = ' Conc max rois'; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

