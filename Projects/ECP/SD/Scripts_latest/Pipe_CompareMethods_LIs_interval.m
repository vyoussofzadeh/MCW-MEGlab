
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

LI_pt = [];
for i=1:length(LI_method_label)
    LI_pt.(LI_method_label{i}) = load(fullfile(data_save_dir,LI_analysis_label{LI_analysis},LI_method_label{i},'LI_Patn'));
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
for i=1:length(LI_method_label)
    LI_pt_val_new.(LI_method_label{i}) = LI_pt.(LI_method_label{i}).LI_sub(:, IA,:);
end

%%
cfg = [];
cfg.thre = fMRI_thre; cfg.LI = fmri_LIs_val;
fmri_LIs_trn = do_ternary_classification2(cfg);
size(fmri_LIs_trn);

%% Plot MEG LI for selected networks
% close all;
network_sel = [1, 2, 6, 11]; % Define the networks to include in the plot
colors = distinguishable_colors(length(network_sel)); % Generate distinct colors for each selected network

figure; % Open a new figure window
hold on; % Keep the plot active to add more elements
plotHandles = gobjects(length(network_sel), 1); % Initialize array for plot handles

% Loop through the selected networks
for net_idx = 1:length(network_sel)
    current_network = network_sel(net_idx); % Current network index

    % Prepare data to plot
    LI_values = []; % Initialize LI_values array
    for i = 1:length(LI_method_label)
        LI_values(i, :) = squeeze(mean(LI_pt_val_new.(LI_method_label{i})(current_network, :, :), 2)); % Extract and average LI values for the current method and network
    end

    % Calculate mean across methods
    meanLI = nanmean(LI_values, 1);
    % Plot the averaged LI values for the current network
    plotHandles(net_idx) = plot(mean(wi'), meanLI, 'LineWidth', 2, 'Color', colors(net_idx,:));
    
    % Find the maximum LI value and its corresponding time
    [maxLI, idx] = max(meanLI);
    maxTime = mean(wi(idx,:));  % Average time at the maximum LI point
    
    % Annotate the maximum value on the plot
    %     text(maxTime, maxLI, sprintf('Mx:%.2f %.2fs', maxLI, maxTime), ...
    %         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    text(maxTime, maxLI, sprintf('%.2fs', maxTime), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    % Draw a vertical line at the max time
    line([maxTime maxTime], ylim, 'Color', colors(net_idx,:), 'LineWidth', 1.5, 'LineStyle', '--');
    
    legendEntries{net_idx} = net_sel_mutiple_label{network_sel(net_idx)}; % Store legend entry
end

xlabel('Time (s)');
ylabel('Laterality Index');
title({'MEG Laterality Index'; 'Over Time for Selected Networks'});
set(gca, 'color', 'none');

% Use the plot handles for the legend to ensure continuity
legend(plotHandles,legendEntries, 'Location', 'southoutside', 'NumColumns', 4, 'Orientation', 'horizontal');

box off;
set(gcf, 'Position', [800, 400, 500, 400]);
hold off; % Release the plot hold

% Set up configuration for exporting the figure
cfg = [];
cfg.outdir = save_dir; % Ensure save_dir is defined and points to a valid directory path
cfg.filename = 'MEG_Laterality_Index_Selected_Networks_Over_Time'; % Filename without the extension
cfg.type = 'svg'; % Specify the type as 'fig'
do_export_fig(cfg); % Call the export function

cd(save_dir); % Change back to the save directory

close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

%% mean MEG li vs. fMRI
% close all,

for i=1:length(LI_method_label)
    
    disp(['correlation analysis: ',LI_method_label{i}])
    
    cfg = []; cfg.wi = wi;
    cfg.ID = sub_MF_pt;
    cfg.thre = MEG_thre;
    cfg.bf = 1;
    cfg.ternary = 0;
    cfg.savefig = 0;
    cfg.outdir = save_dir;
    cfg.net_sel_mutiple_label = net_sel_mutiple_label;
    cfg.net_sel_id = [1,2,6,11];
    cfg.lang_id = {'language_Angular'; 'language_Frontal';'language_Temporal'; 'language_Lateral'};
    cfg.LI_val = LI_pt_val_new.(LI_method_label{i});
    cfg.fmri_LIs_val = fmri_LIs;
    cfg.idx = IB_megfmri;
    cfg.title = LI_method_label{i};
    do_MEG_fMRI_corr_contrast_rois(cfg);
    
    % - export figs
    cfg = []; cfg.outdir = save_dir; filename = ['corr ROIs_', LI_method_label{i}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
    close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

    disp(['Concordance analysis: ',LI_method_label{i}])
    
    disp('---')
    % clc
    cfg = []; cfg.wi = wi;
    cfg.ID = sub_MF_pt;
    cfg.thre = MEG_thre;
    cfg.ternary = 1;
    cfg.savefig = 0;
    cfg.outdir = save_dir;
    cfg.net_sel_mutiple_label = net_sel_mutiple_label;
    cfg.net_sel_id = [1,2,6,11];
    cfg.lang_id = {'language_Angular'; 'language_Frontal';'language_Temporal'; 'language_Lateral'};
    cfg.LI_val = LI_pt_val_new.(LI_method_label{i});
    cfg.fmri_LIs_val = fmri_LIs_trn;
    cfg.idx = IB_megfmri;
    cfg.title = LI_method_label{i};
    cfg.buffervalue = 5;
%     do_MEG_fMRI_concordance_contrast_rois(cfg);
    do_MEG_fMRI_concordance_contrast_rois_interval(cfg);
    
    % - export figs
    cfg = []; cfg.outdir = save_dir; filename = ['concor ROIs_', LI_method_label{i}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
    close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
end
cd(save_dir)

%%
% pause, 

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
    cfg.net_sel_id = [1,2,6,11];
    cfg.lang_id = {'language_Angular'; 'language_Frontal';'language_Temporal'; 'language_Lateral'};
    cfg.LI_val = LI_pt_val_new.(LI_method_label{i});
    cfg.fmri_LIs_val = fmri_LIs;
    cfg.idx = IB_megfmri;
    cfg.title = LI_method_label{i};
    
    % Correlation Analysis
    cfg.thre = MEG_thre;
    cfg.ternary = 0;
    [~, ~, correlationMetrics] = do_MEG_fMRI_corr_contrast_rois2(cfg);
    
    % Concordance Analysis
    cfg.thre = MEG_thre;
    cfg.ternary = 1;
    cfg.buffervalue = 2;
    cfg.fmri_LIs_val = fmri_LIs_trn;
%     [~, ~, concordanceMetrics] = do_MEG_fMRI_concordance_contrast_rois(cfg);
    [~, concordanceMetrics] = do_MEG_fMRI_concordance_contrast_rois_interval(cfg);
    
    % Compile results
    metrics = struct('Correlation', correlationMetrics, 'Concordance', concordanceMetrics);
    resultsTable = [resultsTable; {LI_method_label{i}, metrics}];
end
close all,


%%
% close all,

% Example metric name - adjust according to your actual data structure
metricName = {'Correlation','Concordance'}; % Assuming this is the name of the correlation/concordance field
roi_label = {'Ang', 'Front', 'Temp', 'Lat'};

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

close all

%% Constant 
% pause, 
% close all,

% Assuming metricName contains the names of the fields we want to plot
metricNames = {'Correlation', 'Concordance'};
roi_labels = {'Ang', 'Front', 'Temp', 'Lat'};

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
        maxMethodIndex = 0;

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
                        maxMethodIndex = i;
                        maxTime = mean(wi(localMaxIndex,:));
                    end
                end
            end
        end
        
        if maxValue > -inf
            interval_data = wi(localMaxIndex,:);
            line([maxTime maxTime], ylim, 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--');
            bestMethods{metricIdx, roiIdx} = resultsTable.Method{maxMethodIndex};
            bestMethodName = strtok(bestMethods{metricIdx, roiIdx}, ' ');  %Extracts the first word/name
            text(maxTime + 0.3, maxValue, ...
                sprintf('%s, %.2f; [%.2f,%.2f]', bestMethodName(1), maxValue, interval_data(1), interval_data(end)), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
            %    pause,
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
    set(gcf, 'Position', [1000, 400, 350, 800]);
    
    cfg = []; cfg.outdir = save_dir; filename = [metricNames{metricIdx}, ' rois']; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg);
    close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

end
disp(['saved as, ', filename])

%%
% pause, 
% close all,

% Assuming metricName contains the names of the fields we want to plot
metricNames = {'Correlation', 'Concordance'};
roi_labels = {'Ang', 'Front', 'Temp', 'Lat'};

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

close all

%%
% pause, 
% close all,

% Assuming metricName contains the names of the fields we want to plot
metricNames = {'Correlation', 'Concordance'};
roi_labels = {'Ang', 'Front', 'Temp', 'Lat'};
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
    cfg = []; cfg.outdir = save_dir; filename = [LI_method_labels{methodIdx}, ' Analysis across ROIs']; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
    close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

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
%             maxTime = mean(wi(maxIndex, :), 2); % Assuming 'wi' contains start and end times of intervals
            Timeint = wi(maxIndex, :); % Assuming 'wi' contains start and end times of intervals
            
            % Add to the summary table
            newRow = {LI_method_labels{methodIdx}, metricNames{metricIdx}, roi_labels{roiIdx}, maxValue, Timeint};
            summaryTable = [summaryTable; newRow];
        end
    end
end

% Set column names
summaryTable.Properties.VariableNames = {'LI_Method', 'Metric_Type', 'ROI', 'Max_Value', 'Time_Interval'};

writetable(summaryTable, 'LI_Metrics_Summary.csv');

%
fid = fopen('LI_Metrics_Summary.txt', 'wt'); % Open file for writing
% Print a header
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', summaryTable.Properties.VariableNames{:});

% Loop through each row of the table and print
for i = 1:height(summaryTable)
    fprintf(fid, '%s\t%s\t%s\t%f\t%f\n', summaryTable.LI_Method{i}, summaryTable.Metric_Type{i}, ...
        summaryTable.ROI{i}, summaryTable.Max_Value(i), summaryTable.Time_Interval(i));
end

fclose(fid); % Close the file
disp('constant interval')
disp(summaryTable)


%% Constant interval
% pause, 
% close all,

% Unique ROIs for iteration
uniqueROIs = unique(summaryTable.ROI);

% Loop through each ROI for plotting
figure;
sgtitle('Corr Max'); % Super title for the figure
for i = 1:length(uniqueROIs)
    
    subplot (4,1,i)
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
        text(k, methodData.Max_Value + 0.02, sprintf('%.2f', methodData.Max_Value), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
        k=1+k;
    end
    title([roi]);
    ylabel('Corr.');
    hold off;
    set(gca,'color','none');
    axis tight
    ylim([0, 1])
    set(gcf, 'Position', [100, 100, 200, 800]); % Adjust figure size
end

cfg = []; cfg.outdir = save_dir; filename = ' Corr max'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');


figure;
sgtitle('Conc Max'); % Super title for the figure
for i = 1:length(uniqueROIs)
    
    subplot (4,1, i)
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
        text(k, methodData.Max_Value + 0.02, sprintf('%.2f', methodData.Max_Value), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
        k=1+k;
    end
    title([roi]);
    ylabel('Concordance');
    hold off;
    set(gca,'color','none');
    ylim([0, 100])
    set(gcf, 'Position', [100, 100, 200, 800]); % Adjust figure size    
end

cfg = []; cfg.outdir = save_dir; filename = ' Conc max'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

disp(summaryTable)

%% Dynamic interval analysis
% pause, 
close all,

% Initialize variables
fMRI_LI = fmri_LIs_val;
timePoints = mean(wi, 2);
lowerBound = 0.45; upperBound = 1.1;
lowerBound = 0.4; upperBound = 0.9;
lowerBound = 0.4; upperBound = 0.9;
lowerBound = 0.2; upperBound = 1.1;
lowerBound = 0.35; upperBound = 1.0;
lowerBound = 0.4; upperBound = 1.0;
lowerBound = 0.3; upperBound = 1.2;

% MEG_thre = 0.1;

sub_IDs = sub_MF_pt;
nsub_IDs = cellfun(@(x) [num2str(find(strcmp(sub_IDs, x))), ':', x], sub_IDs, 'UniformOutput', false);

idcx = [1, 2, 6, 11];
MEG_LI_ROIs = cell(1, length(idcx));

for j = 1:length(idcx)
    MEG_LI_ROIs{j} = squeeze(LI_pt_val_new.(LI_method_label{1})(idcx(j), :, :)); % Placeholder for methodIdx
    %     disp(net_sel_mutiple_label{idcx(j)})
end

% fmri_LIs_ROIs = [fmri_LIs.val.language_Angular, fmri_LIs.val.language_Frontal, ...
%     fmri_LIs.val.language_PCingPrecun, fmri_LIs.val.language_Temporal, fmri_LIs.val.language_Lateral];
fmri_LIs_ROIs = [fmri_LIs.val.language_Angular, fmri_LIs.val.language_Frontal, ...
    fmri_LIs.val.language_Temporal, fmri_LIs.val.language_Lateral];
fmri_LIs_ROIs = fmri_LIs_ROIs(IB_megfmri,:);

summaryTableDynamic = table();

for j = 1:length(idcx)
    
    fMRI_LI = fmri_LIs_ROIs(:, j);
    
    for methodIdx = 1:length(LI_method_label)
        
        MEG_LI = squeeze(LI_pt_val_new.(LI_method_label{methodIdx})(idcx(j), :, :));
        
        % Compute group-level correlation
        [groupCorrelation, optimalInterval, optimalInterval_all, pval]= ...
            computeGroupLevelMEGfMRICorrelation_timepoints_interval(MEG_LI, fMRI_LI, wi, lowerBound, upperBound);
        disp(['pval:', num2str(pval)])

%         [groupCorrelation, optimalInterval, ~]= ...
%             computeGroupLevelMEGfMRICorrelation_timepoints_interval(MEG_LI, fMRI_LI, wi, optimalInterval_all(1), optimalInterval_all(2));

        
%         [groupCorrelation, optimalInterval, ~]= computeGroupLevelMEGfMRICorrelation_timepoints_interval(MEG_LI, fMRI_LI, wi, optimalInterval_all(1), optimalInterval_all(2));
%         [groupCorrelation, optimalTimePoints] = computeGroupLevelMEGfMRICorrelation_timepoints(MEG_LI, fMRI_LI, timePoints, lowerBound, upperBound);

        meanOptimalTime = mean(optimalInterval);
        
        [concordance, discordantSubs, groupCorrelation] = ...
            calculateConcordanceForTimePoints_interval(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, wi, optimalInterval);
%         [concordance, discordantSubs] = calculateConcordanceForTimePoints(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimePoints);
        
        % Store results in the summary table
        newRow = {LI_method_labels{methodIdx}, net_sel_mutiple_label{idcx(j)}, groupCorrelation, concordance, meanOptimalTime};
        summaryTableDynamic = [summaryTableDynamic; newRow];
        
        % Plot optimal time points on MEG LI
        % plotOptimalTimePointsOnMEG2(MEG_LI, fMRI_LI, timePoints, optimalTimePoints, discordantSubs, MEG_thre, lowerBound, upperBound);
    end
end

% Set column names for the summary table
summaryTableDynamic.Properties.VariableNames = {'LI_Method', 'ROI', 'Correlation', 'Concordance', 'mean_Optimal_Time'};

% Save summary table
writetable(summaryTableDynamic, 'LI_Metrics_Summary_Dynamic.csv');

% Save dynamic summary table as text file
fid = fopen('LI_Metrics_Summary_Dynamic.txt', 'wt');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', summaryTableDynamic.Properties.VariableNames{:});

for i = 1:height(summaryTableDynamic)
    fprintf(fid, '%s\t%s\t%f\t%f\t%f\n', summaryTableDynamic.LI_Method{i}, summaryTableDynamic.ROI{i}, ...
        summaryTableDynamic.Correlation(i), summaryTableDynamic.Concordance(i), summaryTableDynamic.mean_Optimal_Time(i));
end

fclose(fid);
disp('Summary table for dynamic intervals saved.');
disp(summaryTableDynamic)

%% Dynamic interval plotting
% pause, 
close all,

% Unique ROIs for iteration
uniqueROIs = unique(summaryTableDynamic.ROI);

% Loop through each ROI for plotting
figure;
sgtitle('Corr Max'); % Super title for the figure
for i = 1:length(uniqueROIs)
    
    subplot(4, 1, i)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi), :);
    
    % Plot Max Values for Correlation
    hold on;
    k = 1;
    for method = unique(roiData.LI_Method)'
        methodData = roiData(strcmp(roiData.LI_Method, method), :);
        bar(categorical(methodData.LI_Method), methodData.Correlation, 'BarWidth', 0.2);
        text(k, methodData.Correlation + 0.02, sprintf('%.2f', methodData.Correlation), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
        k = k + 1;
    end
    title(roi);
    ylabel('Corr.');
    hold off;
    set(gca, 'color', 'none');
    axis tight
    ylim([0, 1])
    set(gcf, 'Position', [1000, 400, 200, 700]); % Adjust figure size
end

cfg = []; cfg.outdir = save_dir; filename = 'Corr_dynamic'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');


figure;
sgtitle('Conc Max'); % Super title for the figure
for i = 1:length(uniqueROIs)
    
    subplot(4, 1, i)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi), :);
    
    % Plot Max Values for Concordance
    hold on;
    k = 1;
    for method = unique(roiData.LI_Method)'
        methodData = roiData(strcmp(roiData.LI_Method, method), :);
        bar(categorical(methodData.LI_Method), methodData.Concordance, 'BarWidth', 0.2);
        text(k, methodData.Concordance + 0.02, sprintf('%.2f', methodData.Concordance), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
        k = k + 1;
    end
    title(roi);
    ylabel('Concordance');
    hold off;
    set(gca, 'color', 'none');
    ylim([0, 100])
    set(gcf, 'Position', [1000, 400, 200, 700]); % Adjust figure size
end

cfg = []; cfg.outdir = save_dir; filename = 'Conc_dynamic'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

disp('Dynamic interval analysis plotting completed.');


% Plot difference between Constant and Dynamic Intervals
% pause, 
% close all,

% Extract unique ROIs and methods
uniqueROIs = unique(summaryTable.ROI);
uniqueMethods = unique(summaryTable.LI_Method);

% uniqueROIs_d = unique(summaryTableDynamic.ROI);
% uniqueMethods_d = unique(summaryTableDynamic.LI_Method);

% Initialize difference table
differenceTable = table();

% Compute differences
for i = 1:length(uniqueROIs)
    for j = 1:length(uniqueMethods)
        roi = uniqueROIs{i};
        %         roi_d = uniqueROIs_d{i};
        
        method = uniqueMethods{j};
        %         method_d = uniqueMethods_d{j};
        
        % Extract constant and dynamic metrics
        constantMetrics = summaryTable(strcmp(summaryTable.ROI, roi) & strcmp(summaryTable.LI_Method, method), :);
        dynamicMetrics = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi) & strcmp(summaryTableDynamic.LI_Method, method), :);
        
        if ~isempty(constantMetrics) && ~isempty(dynamicMetrics)
            corrDiff = dynamicMetrics.Correlation - constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Correlation'));
            concDiff = dynamicMetrics.Concordance - constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Concordance'));
            
            % Store in the difference table
            newRow = {method, roi, corrDiff, concDiff};
            differenceTable = [differenceTable; newRow];
        end
    end
end

% Set column names for the difference table
differenceTable.Properties.VariableNames = {'LI_Method', 'ROI', 'Correlation_Diff', 'Concordance_Diff'};

% Plot differences
figure;
sgtitle('Corr Diff, Opt-fixed)');
for i = 1:length(uniqueROIs)
    subplot(4, 1, i)
    roi = uniqueROIs{i};
    roiData = differenceTable(strcmp(differenceTable.ROI, roi), :);
    hold on;
    k = 1;
    for method = unique(roiData.LI_Method)'
        methodData = roiData(strcmp(roiData.LI_Method, method), :);
        bar(categorical(methodData.LI_Method), methodData.Correlation_Diff, 'BarWidth', 0.2);
        text(k, methodData.Correlation_Diff + 0.02, sprintf('%.2f', methodData.Correlation_Diff), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
        k = k + 1;
    end
    title(roi);
    ylabel('Correlation Diff.');
    hold off;
    set(gca, 'color', 'none');
    axis tight
    ylim([-0.5, 0.5])
    set(gcf, 'Position', [1000, 400, 200, 700]); % Adjust figure size
end

cfg = []; cfg.outdir = save_dir; filename = 'Corr_diff'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');


figure;
sgtitle('Conc Diff, Opt-fixed)');
for i = 1:length(uniqueROIs)
    subplot(4, 1, i)
    roi = uniqueROIs{i};
    roiData = differenceTable(strcmp(differenceTable.ROI, roi), :);
    hold on;
    k = 1;
    for method = unique(roiData.LI_Method)'
        methodData = roiData(strcmp(roiData.LI_Method, method), :);
        bar(categorical(methodData.LI_Method), methodData.Concordance_Diff, 'BarWidth', 0.2);
        text(k, methodData.Concordance_Diff + 0.02, sprintf('%.2f', methodData.Concordance_Diff), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
        k = k + 1;
    end
    title(roi);
    ylabel('Concordance Diff.');
    hold off;
    set(gca, 'color', 'none');
    axis tight
    ylim([-30, 30])
    set(gcf, 'Position', [1000, 400, 200, 700]); % Adjust figure size
end

cfg = []; cfg.outdir = save_dir; filename = 'Conc_diff'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');


disp('Difference plotting completed.')

%% Plot Comparison of Constant vs. Dynamic LIs Approach
% close all,

% Extract unique ROIs and methods
uniqueROIs = unique(summaryTable.ROI);
uniqueMethods = unique(summaryTable.LI_Method);

% Initialize comparison table
comparisonTable = table();

% Combine constant and dynamic metrics for comparison
for i = 1:length(uniqueROIs)
    for j = 1:length(uniqueMethods)
        roi = uniqueROIs{i};
        method = uniqueMethods{j};
        
        % Extract constant and dynamic metrics
        constantMetrics = summaryTable(strcmp(summaryTable.ROI, roi) & strcmp(summaryTable.LI_Method, method), :);
        dynamicMetrics = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi) & strcmp(summaryTableDynamic.LI_Method, method), :);
        
        if ~isempty(constantMetrics) && ~isempty(dynamicMetrics)
            % Store in the comparison table
            newRow = {method, roi, 'Constant', constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Correlation')), ...
                constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Concordance'))};
            comparisonTable = [comparisonTable; newRow];
            newRow = {method, roi, 'Dynamic', dynamicMetrics.Correlation, dynamicMetrics.Concordance};
            comparisonTable = [comparisonTable; newRow];
        end
    end
end

% Set column names for the comparison table
comparisonTable.Properties.VariableNames = {'LI_Method', 'ROI', 'Interval_Type', 'Correlation', 'Concordance'};

% Plot comparison
figure;
sgtitle('Correlation: Constant vs. Dynamic Intervals');
for i = 1:length(uniqueROIs)
    subplot(1, length(uniqueROIs), i)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = comparisonTable(strcmp(comparisonTable.ROI, roi), :);
    
    % Prepare data for bar plot
    categories = categorical(roiData.LI_Method);
    intervalTypes = roiData.Interval_Type;
    correlationValues = roiData.Correlation;
    
    % Create bar plot
    barData = reshape(correlationValues, [], 3)';
    b = bar(categories(1:2:end), barData,'BarWidth', 0.2);
    b(1).FaceColor = 'b';
    b(2).FaceColor = 'r';
    
    title(roi);
    ylabel('Correlation');
    set(gca, 'color', 'none');
    ylim([0, 1]);
    box off;
    hold off;
    axis tight
end
lgd = legend(intervalTypes, 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', length(intervalTypes));
lgdPos = lgd.Position; % Get current position
lgdPos(2) = lgdPos(2) - 0.10; % Move legend down
lgd.Position = lgdPos; % Set new position

% legend(intervalTypes, 'Location', 'southoutside', 'Orientation', 'horizontal');
set(gcf, 'Position', [1000, 400, 700, 200]); % Adjust figure size

cfg = []; cfg.outdir = save_dir; filename = 'Correlation_Comparison'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

figure;
sgtitle('Concordance: Constant vs. Dynamic Intervals');
for i = 1:length(uniqueROIs)
    subplot(1, length(uniqueROIs), i)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = comparisonTable(strcmp(comparisonTable.ROI, roi), :);
    
    % Prepare data for bar plot
    categories = categorical(roiData.LI_Method);
    intervalTypes = roiData.Interval_Type;
    concordanceValues = roiData.Concordance;
    
    % Create bar plot
    barData = reshape(concordanceValues, [], 3)';
    b = bar(categories(1:2:end), barData, 'BarWidth', 0.2);
    b(1).FaceColor = 'b';
    b(2).FaceColor = 'r';
    
    title(roi);
    ylabel('Concordance');
    set(gca, 'color', 'none');
    ylim([0, 100]);
    box off;
    hold off;
    axis tight
end
lgd = legend(intervalTypes, 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', length(intervalTypes));
lgdPos = lgd.Position; % Get current position
lgdPos(2) = lgdPos(2) - 0.10; % Move legend down
lgd.Position = lgdPos; % Set new position
set(gcf, 'Position', [100, 400, 700, 200]); % Adjust figure size

cfg = []; cfg.outdir = save_dir; filename = 'Concordance_Comparison'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

disp('Comparison of constant vs. dynamic intervals completed.')


%% Lateral
clc, close all
% Extract unique ROIs and methods
uniqueROIs = unique(summaryTable.ROI);
uniqueMethods = unique(summaryTable.LI_Method);

% Initialize comparison table
comparisonTable = table();

% Combine constant and dynamic metrics for comparison
for i = 3:3%length(uniqueROIs)
    for j = 1:length(uniqueMethods)
        roi = uniqueROIs{i};
        method = uniqueMethods{j};
        
        % Extract constant and dynamic metrics
        constantMetrics = summaryTable(strcmp(summaryTable.ROI, roi) & strcmp(summaryTable.LI_Method, method), :);
        dynamicMetrics = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi) & strcmp(summaryTableDynamic.LI_Method, method), :);
        
        if ~isempty(constantMetrics) && ~isempty(dynamicMetrics)
            % Store in the comparison table
            newRow = {method, roi, 'Constant', constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Correlation')), ...
                constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Concordance'))};
            comparisonTable = [comparisonTable; newRow];
            newRow = {method, roi, 'Dynamic', dynamicMetrics.Correlation, dynamicMetrics.Concordance};
            comparisonTable = [comparisonTable; newRow];
        end
    end
end

% Set column names for the comparison table
comparisonTable.Properties.VariableNames = {'LI_Method', 'ROI', 'Interval_Type', 'Correlation', 'Concordance'};

% Plot comparison
figure;
sgtitle('Corr: Fixed-vs-Opt');
for i = 3:3%length(uniqueROIs)
    subplot(1, 2, 1)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = comparisonTable(strcmp(comparisonTable.ROI, roi), :);
    
    % Prepare data for bar plot
    categories = categorical(roiData.LI_Method);
    intervalTypes = roiData.Interval_Type;
    correlationValues = roiData.Correlation;
    
    % Create bar plot
    barData = reshape(correlationValues, [], 3)';
    b = bar(categories(1:2:end), barData,'BarWidth', 0.2);
    b(1).FaceColor = 'b';
    b(2).FaceColor = 'r';
    disp(barData)
    
    title(roi);
    ylabel('Correlation');
    set(gca, 'color', 'none');
    ylim([0, 1]);
    box off;
    hold off;
    axis tight
end

set(gcf, 'Position', [100, 400, 300, 200]); % Adjust figure size

sgtitle('Conc: Fixed-vs-Opt');
for i = 3:3%length(uniqueROIs)
    subplot(1, 2, 2)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = comparisonTable(strcmp(comparisonTable.ROI, roi), :);
    
    % Prepare data for bar plot
    categories = categorical(roiData.LI_Method);
    intervalTypes = roiData.Interval_Type;
    concordanceValues = roiData.Concordance;
    
    % Create bar plot
    barData = reshape(concordanceValues, [], 3)';
    b = bar(categories(1:2:end), barData, 'BarWidth', 0.2);
    b(1).FaceColor = 'b';
    b(2).FaceColor = 'r';
    disp(barData)
    
    title(roi);
    ylabel('Concordance');
    set(gca, 'color', 'none');
    ylim([0, 100]);
    box off;
    hold off;
    axis tight
end
lgd = legend(intervalTypes, 'Location', 'south', 'Orientation', 'horizontal', 'NumColumns', length(intervalTypes));
set(gcf, 'Position', [400, 400, 300, 250]); % Adjust figure size

cfg = []; cfg.outdir = save_dir; filename = 'Conc_Compr_lat'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

disp('Comparison of constant vs. dynamic intervals completed.');

%% Lateral
clc;
close all;

% Extract unique ROIs and methods
uniqueROIs = unique(summaryTable.ROI);
uniqueMethods = unique(summaryTable.LI_Method);

% Initialize comparison table
comparisonTable = table();

% Combine constant and dynamic metrics for comparison
for i = 1:length(uniqueROIs)
    for j = 1:length(uniqueMethods)
        roi = uniqueROIs{i};
        method = uniqueMethods{j};
        
        % Extract constant and dynamic metrics
        constantMetrics = summaryTable(strcmp(summaryTable.ROI, roi) & strcmp(summaryTable.LI_Method, method), :);
        dynamicMetrics = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi) & strcmp(summaryTableDynamic.LI_Method, method), :);
        
        if ~isempty(constantMetrics) && ~isempty(dynamicMetrics)
            % Store in the comparison table
            newRowConstant = {method, roi, 'Constant', constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Correlation')), ...
                constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Concordance'))};
            newRowDynamic = {method, roi, 'Dynamic', dynamicMetrics.Correlation, dynamicMetrics.Concordance};
            comparisonTable = [comparisonTable; newRowConstant; newRowDynamic];
        end
    end
end

% Set column names for the comparison table
comparisonTable.Properties.VariableNames = {'LI_Method', 'ROI', 'Interval_Type', 'Correlation', 'Concordance'};

% Plot comparison
figure;
sgtitle('Correlation and Concordance: Constant vs. Dynamic Intervals');

for i = 1:length(uniqueROIs)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = comparisonTable(strcmp(comparisonTable.ROI, roi), :);
    
    % Prepare data for bar plot
    categories = categorical(roiData.LI_Method);
    intervalTypes = roiData.Interval_Type;
    correlationValues = roiData.Correlation;
    concordanceValues = roiData.Concordance;
    
    % Create subplots for Correlation and Concordance
    subplot(2, length(uniqueROIs), i);
    barDataCorrelation = reshape(correlationValues, [], 3)';
    bar([1,2], mean(barDataCorrelation), 'grouped');
    title(['Corr - ' roi]);
    ylabel('Correlation');
    ylim([0, 1]);
    set(gca, 'color', 'none');
%     legend(intervalTypes(1:2), 'Location', 'northwest');
    
    subplot(2, length(uniqueROIs), i + length(uniqueROIs));
    barDataConcordance = reshape(concordanceValues, [], 3)';
    bar(1:2, mean(barDataConcordance), 'grouped');
    title(['Con - ' roi]);
    ylabel('Concordance');
    ylim([0, 100]);
    set(gca, 'color', 'none');
end
% legend({'Constant', 'Dynamic'}, 'Location', 'northwest');


% Adjust figure size
set(gcf, 'Position', [1000, 400, 600, 400]);

% Save the figure
cfg = []; 
cfg.outdir = save_dir; 
filename = 'Concordance_Comparison_lat'; 
cfg.filename = filename; 
cfg.type = 'svg'; 
do_export_fig(cfg);

close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

disp('Comparison of constant vs. dynamic intervals completed.');

%% Lateral
clc, close all
% Extract unique ROIs and methods
uniqueROIs = unique(summaryTable.ROI);
uniqueMethods = unique(summaryTable.LI_Method);

% Initialize comparison table
comparisonTable = table();

% Combine constant and dynamic metrics for comparison
for i = 3:3%length(uniqueROIs)
    for j = 1:length(uniqueMethods)
        roi = uniqueROIs{i};
        method = uniqueMethods{j};
        
        % Extract constant and dynamic metrics
        constantMetrics = summaryTable(strcmp(summaryTable.ROI, roi) & strcmp(summaryTable.LI_Method, method), :);
        dynamicMetrics = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi) & strcmp(summaryTableDynamic.LI_Method, method), :);
        
        if ~isempty(constantMetrics) && ~isempty(dynamicMetrics)
            % Store in the comparison table
            newRow = {method, roi, 'Constant', constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Correlation')), ...
                constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Concordance'))};
            comparisonTable = [comparisonTable; newRow];
            newRow = {method, roi, 'Dynamic', dynamicMetrics.Correlation, dynamicMetrics.Concordance};
            comparisonTable = [comparisonTable; newRow];
        end
    end
end

% Set column names for the comparison table
comparisonTable.Properties.VariableNames = {'LI_Method', 'ROI', 'Interval_Type', 'Correlation', 'Concordance'};

% Plot comparison
figure;
sgtitle('Correlation: Constant vs. Dynamic Intervals');
for i = 3:3%length(uniqueROIs)
    subplot(1, 2, 1)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = comparisonTable(strcmp(comparisonTable.ROI, roi), :);
    
    % Prepare data for bar plot
    categories = categorical(roiData.LI_Method);
    intervalTypes = roiData.Interval_Type;
    correlationValues = roiData.Correlation;
    
    % Create bar plot
    barData = reshape(correlationValues, [], 3)';
    b = bar(mean(barData,1),'BarWidth', 0.2);
    b(1).FaceColor = 'b';
%     b(2).FaceColor = 'r';
    
    title(roi);
    ylabel('Correlation');
    set(gca, 'color', 'none');
    ylim([0, 1]);
    box off;
    hold off;
    axis tight
end

set(gcf, 'Position', [100, 400, 300, 200]); % Adjust figure size

sgtitle('Concordance: Constant vs. Dynamic Intervals');
for i = 3:3%length(uniqueROIs)
    subplot(1, 2, 2)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = comparisonTable(strcmp(comparisonTable.ROI, roi), :);
    
    % Prepare data for bar plot
    categories = categorical(roiData.LI_Method);
    intervalTypes = roiData.Interval_Type;
    concordanceValues = roiData.Concordance;
    
    % Create bar plot
    barData = reshape(concordanceValues, [], 3)';
    b = bar(mean(barData,1), 'BarWidth', 0.2);
    b(1,1).FaceColor = 'b';
%     b(2,1).FaceColor = 'r';
    
    title(roi);
    ylabel('Concordance');
    set(gca, 'color', 'none');
    ylim([0, 100]);
    box off;
    hold off;
    axis tight
end
lgd = legend(intervalTypes, 'Location', 'south', 'Orientation', 'horizontal', 'NumColumns', length(intervalTypes));
set(gcf, 'Position', [400, 400, 300, 250]); % Adjust figure size

cfg = []; cfg.outdir = save_dir; filename = 'Concordance_Comparison_lat'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

disp('Comparison of constant vs. dynamic intervals completed.')
