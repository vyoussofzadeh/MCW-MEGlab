
%
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
sgtitle('Corr Diff, Opt-fixed');
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
sgtitle('Conc Diff, Opt-fixed');
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