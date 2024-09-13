function plotSummaryAndDifferences(summaryTableDynamic, save_dir)
    % Define unique ROIs and methods
    uniqueROIs = unique(summaryTableDynamic.ROI);
    uniqueMethods = unique(summaryTableDynamic.LI_Method);
    
    % Initialize a difference table for storing metric differences
    differenceTable = table();

    % Open a new figure
    figure;
    set(gcf, 'Position', [1000, 400, 400, 1000]); % Larger figure to accommodate all subplots

    % Iterate over each ROI
    for i = 1:length(uniqueROIs)
        roi = uniqueROIs{i};
        
        % Extract data for the current ROI
        roiData = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi), :);

        % Plot Maximum Correlation values
        subplot(8, 1, 2*i-1);
        hold on;
        k = 1;
        for method = unique(roiData.LI_Method)'
            methodData = roiData(strcmp(roiData.LI_Method, method), :);
            bar(categorical(methodData.LI_Method), methodData.Correlation, 'BarWidth', 0.2);
            text(k, methodData.Correlation + 0.02, sprintf('%.2f', methodData.Correlation), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
            k = k + 1;
        end
        
        
        
%         methodData = roiData(strcmp(roiData.LI_Method, method), :);
%         bar(categorical(roiData.LI_Method), roiData.Correlation, 'BarWidth', 0.2);
%         text(i, roiData.Correlation + 0.02, sprintf('%.2f', roiData.Correlation), ...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
        ylabel('Max Corr.');
        title(['Corr Max - ', roi]);
        set(gca, 'color', 'none');
        ylim([0, 1]);

        % Plot Maximum Concordance values
        subplot(8, 1, 2*i);
        hold on;
        k = 1;
        for method = unique(roiData.LI_Method)'
            methodData = roiData(strcmp(roiData.LI_Method, method), :);
            bar(categorical(methodData.LI_Method), methodData.Concordance, 'BarWidth', 0.2);
            text(k, methodData.Correlation + 0.02, sprintf('%.2f', methodData.Concordance), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
            k = k + 1;
        end
        
%         bar(categorical(roiData.LI_Method), roiData.Concordance, 'BarWidth', 0.2);
        ylabel('Max Conc.');
        title(['Conc Max - ', roi]);
        set(gca, 'color', 'none');
        ylim([0, 100]);

        % Calculate differences and add to the difference table
        for j = 1:length(uniqueMethods)
            method = uniqueMethods{j};
            constantMetrics = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi) & strcmp(summaryTableDynamic.LI_Method, method), :);
            dynamicMetrics = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi) & strcmp(summaryTableDynamic.LI_Method, method), :);

            if ~isempty(constantMetrics) && ~isempty(dynamicMetrics)
                corrDiff = dynamicMetrics.Correlation - constantMetrics.Correlation;
                concDiff = dynamicMetrics.Concordance - constantMetrics.Concordance;
                differenceTable = [differenceTable; {method, roi, corrDiff, concDiff}];
            end
        end
    end

    % Label the figure
    sgtitle('Summary Metrics - Dynamic Intervals Analysis');
    
    % Save the figure
    saveas(gcf, fullfile(save_dir, 'Combined_Metrics_Summary.fig'));
    disp('Combined metrics plot saved.');

    % Optionally, print the difference table to the command window
    disp(differenceTable);
end
