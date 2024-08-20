%% Lateral
% clc;
% close all;

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
    box off
    set(gca, 'color', 'none');
%     legend(intervalTypes(1:2), 'Location', 'northwest');
    
    subplot(2, length(uniqueROIs), i + length(uniqueROIs));
    barDataConcordance = reshape(concordanceValues, [], 3)';
    bar(1:2, mean(barDataConcordance), 'grouped');
    box off
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
filename = 'Concordance_Comparison_lat1'; 
cfg.filename = filename; 
cfg.type = 'svg'; 
do_export_fig(cfg);

close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

disp('Comparison of constant vs. dynamic intervals completed.');
