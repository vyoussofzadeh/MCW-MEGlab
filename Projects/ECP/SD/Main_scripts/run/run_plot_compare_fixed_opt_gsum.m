%% Clear Workspace and Close Figures
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
            newRow = {method, roi, 'Constant', mean(constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Correlation'))), ...
                mean(constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Concordance')))};
            comparisonTable = [comparisonTable; newRow];
            newRow = {method, roi, 'Dynamic', mean(dynamicMetrics.Correlation), mean(dynamicMetrics.Concordance)};
            comparisonTable = [comparisonTable; newRow];
        end
    end
end

% Set column names for the comparison table
comparisonTable.Properties.VariableNames = {'LI_Method', 'ROI', 'Interval_Type', 'Correlation', 'Concordance'};

% Plot comparison
figure;
sgtitle('Comparison: Constant vs. Dynamic Intervals');

% Correlation comparison
subplot(1, 2, 1)
barData = zeros(length(uniqueMethods), 2);
for j = 1:length(uniqueMethods)
    method = uniqueMethods{j};
    methodData = comparisonTable(strcmp(comparisonTable.LI_Method, method), :);
    barData(j, :) = [mean(methodData.Correlation(strcmp(methodData.Interval_Type, 'Constant'))), ...
                     mean(methodData.Correlation(strcmp(methodData.Interval_Type, 'Dynamic')))];
end
b = bar(categorical(uniqueMethods), barData, 'BarWidth', 0.4);
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
title('Correlation');
ylabel('Value');
ylim([0, 1]);
box off;
set(gca,'color','none');

legend('Constant', 'Dynamic', 'Location', 'southoutside', 'Orientation', 'horizontal');

% Concordance comparison
subplot(1, 2, 2)
barData = zeros(length(uniqueMethods), 2);
for j = 1:length(uniqueMethods)
    method = uniqueMethods{j};
    methodData = comparisonTable(strcmp(comparisonTable.LI_Method, method), :);
    barData(j, :) = [mean(methodData.Concordance(strcmp(methodData.Interval_Type, 'Constant'))), ...
                     mean(methodData.Concordance(strcmp(methodData.Interval_Type, 'Dynamic')))];
end
b = bar(categorical(uniqueMethods), barData, 'BarWidth', 0.4);
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
title('Concordance');
ylabel('Value');
ylim([0, 100]);
box off;
set(gca,'color','none');


legend('Constant', 'Dynamic', 'Location', 'southoutside', 'Orientation', 'horizontal');
set(gcf, 'Position', [100, 100, 700, 300]); % Adjust figure size

% Export figure
cfg = struct('outdir', save_dir, 'filename', 'Comparison_Constant_Dynamic', 'type', 'svg');
do_export_fig(cfg);
close all;
combined_path = fullfile(save_dir, [cfg.filename, '.svg']);
web(combined_path, '-new');

disp('Comparison of constant vs. dynamic intervals completed.');

% %% Lateral
% % clc, close all
% % Extract unique ROIs and methods
% uniqueROIs = unique(summaryTable.ROI);
% uniqueMethods = unique(summaryTable.LI_Method);
% 
% % Initialize comparison table
% comparisonTable = table();
% 
% % Combine constant and dynamic metrics for comparison
% for i = 3:3%length(uniqueROIs)
%     for j = 1:length(uniqueMethods)
%         roi = uniqueROIs{i};
%         method = uniqueMethods{j};
%         
%         % Extract constant and dynamic metrics
%         constantMetrics = summaryTable(strcmp(summaryTable.ROI, roi) & strcmp(summaryTable.LI_Method, method), :);
%         dynamicMetrics = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi) & strcmp(summaryTableDynamic.LI_Method, method), :);
%         
%         if ~isempty(constantMetrics) && ~isempty(dynamicMetrics)
%             % Store in the comparison table
%             newRow = {method, roi, 'Constant', constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Correlation')), ...
%                 constantMetrics.Max_Value(strcmp(constantMetrics.Metric_Type, 'Concordance'))};
%             comparisonTable = [comparisonTable; newRow];
%             newRow = {method, roi, 'Dynamic', dynamicMetrics.Correlation, dynamicMetrics.Concordance};
%             comparisonTable = [comparisonTable; newRow];
%         end
%     end
% end
% 
% % Set column names for the comparison table
% comparisonTable.Properties.VariableNames = {'LI_Method', 'ROI', 'Interval_Type', 'Correlation', 'Concordance'};
% 
% % Plot comparison
% figure;
% sgtitle('Correlation: Constant vs. Dynamic Intervals');
% for i = 3:3%length(uniqueROIs)
%     subplot(1, 2, 1)
%     roi = uniqueROIs{i};
%     
%     % Extract data for the current ROI
%     roiData = comparisonTable(strcmp(comparisonTable.ROI, roi), :);
%     
%     % Prepare data for bar plot
%     categories = categorical(roiData.LI_Method);
%     intervalTypes = roiData.Interval_Type;
%     correlationValues = roiData.Correlation;
%     
%     % Create bar plot
%     barData = reshape(correlationValues, [], 3)';
%     b = bar(mean(barData,1),'BarWidth', 0.2);
% %     b(1).FaceColor = 'b';
% %     b(2).FaceColor = 'r';
%     
%     title(roi);
%     ylabel('Correlation');
%     set(gca, 'color', 'none');
%     ylim([0, 1]);
%     box off;
%     hold off;
%     axis tight
% end
% 
% set(gcf, 'Position', [100, 400, 300, 200]); % Adjust figure size
% 
% sgtitle('Concordance: Constant vs. Dynamic Intervals');
% for i = 3:3%length(uniqueROIs)
%     subplot(1, 2, 2)
%     roi = uniqueROIs{i};
%     
%     % Extract data for the current ROI
%     roiData = comparisonTable(strcmp(comparisonTable.ROI, roi), :);
%     
%     % Prepare data for bar plot
%     categories = categaorical(roiData.LI_Method);
%     intervalTypes = roiData.Interval_Type;
%     concordanceValues = roiData.Concordance;
%     
%     % Create bar plot
%     barData = reshape(concordanceValues, [], 3)';
%     b = bar(mean(barData,1), 'BarWidth', 0.2);
% %     b(1,1).FaceColor = 'b';
% %     b(2,1).FaceColor = 'r';
%     
%     title(roi);
%     ylabel('Concordance');
%     set(gca, 'color', 'none');
%     ylim([0, 100]);
%     box off;
%     hold off;
%     axis tight
% end
% 
% lgd = legend(intervalTypes, 'Location', 'southout', 'Orientation', 'horizontal', 'NumColumns', length(intervalTypes));
% set(gcf, 'Position', [400, 400, 300, 250]); % Adjust figure size
% 
% cfg = []; cfg.outdir = save_dir; filename = 'Concordance_Comparison_lat2'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
% close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
% 
% disp('Comparison of constant vs. dynamic intervals completed.')
