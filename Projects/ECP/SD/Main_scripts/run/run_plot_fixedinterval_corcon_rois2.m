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
    cfg = []; cfg.outdir = save_dir; filename = [[roi_labels{roiIdx}], ' CorConc cmpr']; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
    close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
end

% close all