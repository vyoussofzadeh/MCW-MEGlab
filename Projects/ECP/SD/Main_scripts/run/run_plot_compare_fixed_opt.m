% clc, close all
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
    b = bar(categories(1:2:end), barData,'BarWidth', 0.5);
    b(1).FaceColor = [.6 .6 .6];%'b';
    b(2).FaceColor = rgb(66, 165, 245) ; %'r';
    disp(barData)
    
    title(roi);
    ylabel('Correlation');
    set(gca, 'color', 'none');
    ylim([0, 1]);
    box off;
    hold off;
    %     axis tight
end

lgd = legend(intervalTypes, 'Location', 'bestoutside', 'Orientation', 'horizontal', 'NumColumns', length(intervalTypes));
lgd.Visible = 'off';   % Hides the legend but does not delete it

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
    b = bar(categories(1:2:end), barData, 'BarWidth', 0.5);
    b(1).FaceColor = [.6 .6 .6];%'b';
    b(2).FaceColor = rgb(66, 165, 245) ; %'r';
    disp(barData)
    ylim([0, 100]);
    title(roi);
    ylabel('LI Concordance (%)');
    set(gca, 'color', 'none');
    box off;
    hold off;
    %     axis tight
end

if metricIdx == 2
    lgd.Visible = 'off';   % Hides the legend but does not delete it
end

lgd = legend(intervalTypes, 'Location', 'bestoutside', 'Orientation', 'horizontal', 'NumColumns', length(intervalTypes));

% lgd = legend(intervalTypes, 'Location', 'south', 'Orientation', 'horizontal', 'NumColumns', length(intervalTypes));
set(gcf, 'Position', [400, 400, 300, 350]); % Adjust figure size

cfg = []; cfg.outdir = save_dir; filename = 'Conc_Compr_lat'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

disp('Comparison of constant vs. dynamic intervals completed.');