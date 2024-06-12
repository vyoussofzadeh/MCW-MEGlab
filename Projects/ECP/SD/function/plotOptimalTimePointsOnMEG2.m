function plotOptimalTimePointsOnMEG2(MEG_LI, fMRI_LI, timePoints, optimalTimePoints, discordantSamples, MEG_thre, lowerBound, upperBound)

numSubjects = size(MEG_LI, 1);

% Create a figure to hold all subplots
figure;

for subj = 1:numSubjects
    % Find the index of the optimal time point
    optimalTimePointIdx = find(timePoints == optimalTimePoints(subj));
    
    subplot(ceil(numSubjects / 6), 6, subj); % Adjust for 6 columns
    plot(timePoints, MEG_LI(subj, :), 'b'); % Plot MEG LI as a blue line
    hold on;
    
    % Plot grey lines for the boundaries
    yline(100*MEG_thre, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); % Grey line for upper boundary
    yline(-100*MEG_thre, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); % Grey line for lower boundary
    
    % Plot grey lines for the boundaries
    xline(lowerBound, 'Color', [0.5 0 0.5], 'LineStyle', '-'); % Grey line for upper boundary
    xline(upperBound, 'Color', [0.5 0 0.5], 'LineStyle', '-'); % Grey line for lower boundary
    
    % Check if the current subject is a discordant sample
    if ismember(subj, discordantSamples)
        % Highlight the optimal time point for discordant samples
        plot(timePoints(optimalTimePointIdx), MEG_LI(subj, optimalTimePointIdx), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); % Red circle for discordant samples
    else
        % Highlight the optimal time point for non-discordant samples
        plot(timePoints(optimalTimePointIdx), MEG_LI(subj, optimalTimePointIdx), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'g'); % Green circle for non-discordant samples
    end
    
    % Set the y-axis limits
    ylim([-100 100]);
    
    hold off;
    set(gca, 'color', 'none');
    title(sprintf('S%d|fMRI:%.1f|MEG:%.1f', subj, fMRI_LI(subj), MEG_LI(subj, optimalTimePointIdx)));
end

xlabel('Time Points');
ylabel('MEG LI');
set(gcf, 'Position', [100 300 1200 1000]);

end

