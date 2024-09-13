function plotOptimalTimePointsOnMEG3(rSNR_MEG, MEG_LI, fMRI_LI, wi, optimalIndices, discordantSamples, MEG_thre, lowerBound, upperBound)

numSubjects = size(MEG_LI, 1);

% timePoints = mean(wi,2);
timePoints = wi(:,1);

% Create a figure to hold all subplots
figure;

for subj = 1:numSubjects
    % Find the index of the optimal time point
    optimalTimePointIdx = optimalIndices(subj);
    
    subplot(ceil(numSubjects / 6), 6, subj); % Adjust for 6 columns
    yyaxis left; % Left Y-axis for original SNR values
    p1 = plot(timePoints, MEG_LI(subj, :), 'DisplayName', 'MEG LI', 'LineWidth', 1.5); % Plot MEG LI as a blue line
    hold on;
    
    % Plot grey lines for the boundaries
    yline(100*MEG_thre, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); % Grey line for upper boundary
    yline(-100*MEG_thre, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); % Grey line for lower boundary
    
    % Plot grey lines for the boundaries
    xline(lowerBound, 'Color', [0.5 0 0.5], 'LineStyle', '-'); % Grey line for upper boundary
    xline(upperBound, 'Color', [0.5 0 0.5], 'LineStyle', '-'); % Grey line for lower boundary
    
    xlim([timePoints(1) timePoints(end)])
%     % Check if the current subject is a discordant sample
%     if ismember(subj, discordantSamples)
%         % Highlight the optimal time point for discordant samples
%         plot(timePoints(optimalTimePointIdx), MEG_LI(subj, optimalTimePointIdx), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); % Red circle for discordant samples
%     else
%         % Highlight the optimal time point for non-discordant samples
%         plot(timePoints(optimalTimePointIdx), MEG_LI(subj, optimalTimePointIdx), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'g'); % Green circle for non-discordant samples
%     end
%     
    % Set the y-axis limits
    ylim([-100 100]);
    box off
    hold on;
    set(gca, 'color', 'none');
    title(sprintf('S%d|fMRI:%.1f|MEG:%.1f', subj, fMRI_LI(subj), MEG_LI(subj, optimalTimePointIdx)));
    
    
    currentDiff = abs((rSNR_MEG.rSNR_left(subj, :)) - (rSNR_MEG.rSNR_right(subj, :)));
    yyaxis right; % Left Y-axis for original SNR values
    p2 = plot(timePoints, currentDiff, 'DisplayName', 'rSNR diff', 'LineWidth', 1.0, 'LineStyle', '--'); % Plot left SNR in green
    
    % Check if the current subject is a discordant sample
    if ismember(subj, discordantSamples)
        % Highlight the optimal time point for discordant samples
        p3 = plot(timePoints(optimalTimePointIdx), currentDiff(optimalTimePointIdx), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Optimal time point'); % Red circle for discordant samples
    else
        % Highlight the optimal time point for non-discordant samples
        p3 = plot(timePoints(optimalTimePointIdx), currentDiff(optimalTimePointIdx), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'DisplayName', 'Optimal time point'); % Green circle for non-discordant samples
    end
    
    
end

yyaxis left; % Switch back to left axis for common settings
xlabel('Time Points');
ylabel('MEG LI');

yyaxis right; % Switch back to left axis for common settings
ylabel('rSNR abs(L-R)');

set(gcf, 'Position', [100, 100, 1600, 1300]);

lgd = legend([p1, p2, p3], 'Location', 'best','Orientation', 'horizontal'); % Add a legend
lgdPos = lgd.Position; % Get current position
lgdPos(2) = lgdPos(2) - 0.09; % Move legend down
lgdPos(1) = lgdPos(1) - 0.05; % Move legend down
lgd.Position = lgdPos;


end

