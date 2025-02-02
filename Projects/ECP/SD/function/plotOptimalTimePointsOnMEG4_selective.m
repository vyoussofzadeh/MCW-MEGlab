function plotOptimalTimePointsOnMEG4_selective(rSNR_MEG, MEG_LI, fMRI_LI, wi, optimalIndices, discordantSamples, MEG_thre, bounds)

numSubjects = length(discordantSamples);

% timePoints = mean(wi,2);
timePoints = wi(:,1);

% Create a figure to hold all subplots
figure;

for subj = 1:length(discordantSamples)
    % Find the index of the optimal time point
    
    subj_sel = discordantSamples(subj);
    
    optimalTimePointIdx = optimalIndices(subj_sel);
    
    %     subplot(ceil(numSubjects / 6), 6, subj); % Adjust for 6 columns
    
    if length(discordantSamples) > 10
        subplot(ceil(numSubjects / 6), 6, subj); % Adjust for 6 columns
        set(gcf, 'Position', [100, 100, 1600, 1300]);
    else
        subplot(ceil(numSubjects / 2), 2, subj); % Adjust for 6 columns
        set(gcf, 'Position', [200, 400, 800, 500]);
    end
    
    yyaxis left; % Left Y-axis for original SNR values
    p1 = plot(timePoints, MEG_LI(subj_sel, :), 'DisplayName', 'MEG LI','LineWidth', 1.5); % Plot MEG LI as a blue line
    hold on;
    
    % Plot grey lines for the boundaries
    yline(100*MEG_thre, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); % Grey line for upper boundary
    yline(-100*MEG_thre, 'Color', [0.5 0.5 0.5], 'LineStyle', '--'); % Grey line for lower boundary
    
    % Plot grey lines for the boundaries
    xline(timePoints(bounds(subj,1)), 'Color', [0.5 0 0.5], 'LineStyle', '-'); % Grey line for upper boundary
    xline(timePoints(bounds(subj,2)), 'Color', [0.5 0 0.5], 'LineStyle', '-'); % Grey line for lower boundary
    
    % Set the y-axis limits
    ylim([-100 100]);
    box off
%     hold on;
    set(gca, 'color', 'none');
    title(sprintf('S%d|fMRI:%.1f|MEG:%.1f', subj_sel, fMRI_LI(subj_sel), MEG_LI(subj_sel, optimalTimePointIdx)));
    
    
    currentDiff = rSNR_MEG.rSNR_left(subj_sel, :) - rSNR_MEG.rSNR_right(subj_sel, :);
    yyaxis right; % Left Y-axis for original SNR values
    p2 = plot(timePoints, currentDiff, 'DisplayName', 'rSNR diff', 'LineWidth', 1.5, 'LineStyle', '--'); % Plot left SNR in green
    
    % Check if the current subject is a discordant sample
    if ismember(subj_sel, discordantSamples)
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
ylabel('rSNR (L - R)');
% set(gcf, 'Position', [100 300 1200 1000]);

% hold off
lgd = legend([p1, p2, p3], 'Location', 'best','Orientation', 'horizontal'); % Add a legend
lgdPos = lgd.Position; % Get current position
lgdPos(2) = lgdPos(2) - 0.22; % Move legend down
lgd.Position = lgdPos;

end

