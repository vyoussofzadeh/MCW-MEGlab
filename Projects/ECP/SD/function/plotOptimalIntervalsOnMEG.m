function plotOptimalIntervalsOnMEG(MEG_LI, fMRI_LI, ID, timePoints, optimalIntervals)

    numSubjects = size(MEG_LI, 1);
    
    % Create a figure to hold all subplots
    figure;
    
    for subj = 1:numSubjects
        % Calculate the mean MEG LI within the optimal interval
        intervalStartIdx = find(timePoints == optimalIntervals(subj, 1));
        intervalEndIdx = find(timePoints == optimalIntervals(subj, 2));
        meanOptimalMEG_LI = mean(MEG_LI(subj, intervalStartIdx:intervalEndIdx));
        
        subplot(ceil(numSubjects / 5), 5, subj); % Adjust for 5 columns
        plot(timePoints, MEG_LI(subj, :), 'b'); % Plot MEG LI as a blue line
        hold on;
        
        % Highlight the optimal interval
        patch([timePoints(intervalStartIdx), timePoints(intervalEndIdx), timePoints(intervalEndIdx), timePoints(intervalStartIdx)], ...
              [min(MEG_LI(subj, :)), min(MEG_LI(subj, :)), max(MEG_LI(subj, :)), max(MEG_LI(subj, :))], ...
              'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Adjust color and transparency as needed

        % Set the y-axis limits
        ylim([-100 100]);

        hold off;
        title(sprintf('Subject %s | fMRI LI: %.2f | Mean MEG LI: %.2f', ID{subj}, fMRI_LI(subj), meanOptimalMEG_LI));
        xlabel('Time Points');
        ylabel('MEG LI');
    end
end

