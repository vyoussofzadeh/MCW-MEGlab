% Assuming metricName contains the names of the fields we want to plot
metricNames = {'Correlation', 'Concordance'};
roi_labels = {'Ang', 'Front', 'Temp', 'Lat'};
% LI_method_labels = {'SourceMag', 'Count', 'Bootstrp'}; % LI Methods
LI_method_labels = LI_method_label;

% Colors for each ROI, adjust or extend as needed
colors = lines(length(roi_labels));

colors = [
    0.96 0.49 0
    0.22 0.56 0.24
    0.69 0.71 0.17
    0.48 0.12 0.66
    ];

colors = [
    0 114 189
    217 83 25
    237 177 32
    126 47 142
    ]/256;

midpoints = mean(wi, 2);

% Loop through each LI method
for methodIdx = 1:length(LI_method_labels)
    figure; % Initialize a figure for the current LI method
    %     figure('Color','w','Units','normalized','Position',[0.4 0.2 0.25 0.6]);
    
    sgtitle([LI_method_labels{methodIdx}, ' analysis: ROIs']); % Super title for the figure
    %     set(gcf, 'Position', [100, 100, 600, 600]); % Adjust figure size
    
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
                    plot(midpoints, currentMetrics(roiIdx, :), 'LineWidth', 3, 'Color', colors(roiIdx,:));
                    
                    % Set x-ticks and labels for every other interval
%                     xticks(midpoints(1:3:end));  % x-ticks for every other midpoint
%                     xLabels = arrayfun(@(i) sprintf('%.1f - %.1f', wi(i, 1), wi(i, 2)), 1:length(midpoints), 'UniformOutput', false);
%                     xticklabels(xLabels(1:3:end));  % x-labels for every other interval
%                     xtickangle(45);  % Rotate labels by 45 degrees
                    
                    % Set font size for x-axis labels
                    ax = gca; % Get current axes
                    ax.XAxis.FontSize = 8; % Set font size for x-axis only
                    ax.YAxis.FontSize = 8; % Set font size for y-axis only
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
        
        if metricIdx == 2
            lgd.Visible = 'off';   % Hides the legend but does not delete it
        end
        
        lgd = legend(roi_labels, 'Location', 'bestoutside', 'Orientation', 'horizontal', 'NumColumns', length(roi_labels));
        %         lgdPos = lgd.Position; % Get current position
        %         lgdPos(2) = lgdPos(2) - 0.1; % Move legend down
        %         lgd.Position = lgdPos; % Set new position
        
        
        hold off; % Release hold for next metric type subplot
    end
    cfg = []; cfg.outdir = save_dir; filename = [LI_method_labels{methodIdx}, ' Analysis across ROIs']; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
    close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    
    
end