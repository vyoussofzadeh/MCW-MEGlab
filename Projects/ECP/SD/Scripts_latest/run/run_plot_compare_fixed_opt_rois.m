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
