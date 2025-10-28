% Assuming metricName contains the names of the fields we want to plot
metricNames = {'Correlation', 'Concordance'};
roi_labels = {'Ang', 'Front', 'Temp', 'Lat'};

baseColors = [
    31, 78, 121;
    132, 60, 12;
    127, 96, 0;
    112, 48, 160;
    84, 130, 53]/ 255;

% baseColors = [
%     0.96 0.49 0
%     0.22 0.56 0.24
%     0.69 0.71 0.17
%     0.48 0.12 0.66
%     ];

methodColors = [
    0.00 0.45 0.74;   % Method 1  (e.g., Source-Magnitude)
    0.85 0.33 0.10;   % Method 2  (e.g., Counting-Vertices)
    0.47 0.67 0.19];  % Method 3  (e.g., Bootstrapping)

methodColors = [ ...
    0.17 0.35 0.73;   % Method 1  indigo-blue
    0.87 0.48 0.09;   % Method 2  amber-orange
    0.00 0.60 0.45];  % Method 3  jade-green

methodColors = [ ...
    0.96 0.49 0;   % Method 1  indigo-blue
    0.69 0.71 0.17;   % Method 2  amber-orange
    0.48 0.12 0.66];  % Method 3  jade-green

methodColors = lines(3);              %  [0 0.4470 0.7410;

methodColors = [ ...
    0.60 0.60 0.60;    % Method-1  ? grey
    0.83 0.25 0.25;    % Method-2  ? red   (#D33F3F-ish)
    0.15 0.44 0.76];   % Method-3  ? blue  (#276FD1-ish)

% Enhance color contrast
contrastFactorMagnitude = 0.5;  % Higher factor for lighter colors in Magnitude
contrastFactorCounting = 0.2;   % Lower factor for slightly brighter colors in Counting

colorsMagnitude = min(baseColors + contrastFactorMagnitude, 1);
colorsCounting = min(baseColors + contrastFactorCounting, 1);

midpoints = mean(wi, 2);


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
            
%             if i > size(methodColors,1)
%                 colorToUse = [0.6 0.6 0.6];   % fallback if you have >3 methods
%             else
%                 colorToUse = methodColors(i,:);
%             end
            
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
                    
                    plot(midpoints, dataToPlot, 'LineWidth', 2, 'Color', colorToUse);
                    
                    
                    % Set x-ticks and labels for every other interval
                    xticks(midpoints(1:end));  % x-ticks for every other midpoint
                    %                     xLabels = arrayfun(@(i) sprintf('%.1f - %.1f', wi(i, 1), wi(i, 2)), 1:length(midpoints), 'UniformOutput', false);
                    %                     xticklabels(xLabels(1:3:end));  % x-labels for every other interval
                    xtickangle(90);  % Rotate labels by 45 degrees
                    axis tight
                    
                    % Set font size for x-axis labels
                    ax = gca; % Get current axes
                    ax.XAxis.FontSize = 8; % Set font size for x-axis only
                    ax.YAxis.FontSize = 8; % Set font size for y-axis only
                    
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
            line([maxTime maxTime], ylim, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');
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
            ylim([50 90]);
        end
        
        if roiIdx == length(roi_labels)
            %             ylabel(['LIs ', metricNames{metricIdx}, ' (MEG vs. fMRI)']);
            ylabel([metricNames{metricIdx}]);
            xlabel('Time (sec)');
            %         end
        else
            % every other subplot ? no x-tick labels & no legend
            set(gca,'XTickLabel',[]);   % keeps ticks but removes numbers
            %  (use `set(gca,'XTick',[])` if you want the ticks gone too)
        end
        lgd = legend(resultsTable.Method, 'Location', 'southout', 'NumColumns', 2, 'Orientation', 'horizontal');
        
        title(roi_labels{roiIdx});
        box off;
        %         axis square
        hold off;
    end
    
    % Super title for the entire figure
    sgtitle([metricNames{metricIdx}, ' Analysis']);
    set(gcf, 'Position', [1000, 400, 350, 1100]);
    
    cfg = []; cfg.outdir = save_dir; filename = [metricNames{metricIdx}, ' rois']; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg);
    close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    
end
disp(['saved as, ', filename])

