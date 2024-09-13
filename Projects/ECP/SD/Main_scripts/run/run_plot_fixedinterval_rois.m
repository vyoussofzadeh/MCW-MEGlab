% Assuming metricName contains the names of the fields we want to plot
metricNames = {'Correlation', 'Concordance'};
roi_labels = {'Ang', 'Front', 'Temp', 'Lat'};
% LI_method_labels = {'SourceMag', 'Count', 'Bootstrp'}; % LI Methods
LI_method_labels = LI_method_label;

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