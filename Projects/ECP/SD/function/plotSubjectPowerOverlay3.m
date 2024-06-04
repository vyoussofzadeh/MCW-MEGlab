function [optimalTimePoints, interval] = plotSubjectPowerOverlay3(powerLeft, powerRight, ID, timePoints, plot_flag, lowerBound, upperBound)
% Plots overlaid power values for left and right hemispheres for all subjects with individual peak annotations.
% Searches for the optimal time point based on the combined power values of the left and right hemispheres within the interval defined by lowerBound and upperBound.
% Additionally, plots the distribution of peaks for left and right power source values and the distribution of latency peaks across subjects if plotting is enabled.

numSubjects = size(powerLeft, 1);
optimalTimePoints = zeros(numSubjects, 1); % Initialize array to store optimal time points

% Create a figure to hold all subplots if plotting is enabled
if plot_flag, figure; end

for subj = 1:numSubjects
    combinedPower = (powerLeft(subj, :) - powerRight(subj, :))./(powerLeft(subj, :) + powerRight(subj, :)); % Calculate combined power
    combinedPower = abs(combinedPower);
    optimalPower = -inf; % Initialize maximum power
    optimalIdx = NaN; % Initialize index of optimal time point
    
    % Loop over all time points within the interval
    for i = 1:length(timePoints)
        if timePoints(i) >= lowerBound && timePoints(i) <= upperBound && combinedPower(i) > optimalPower
            optimalPower = combinedPower(i);
            optimalIdx = i;
        end
    end
    
    % Store the optimal time point
    optimalTimePoints(subj) = timePoints(optimalIdx);
    
    if plot_flag
        % Plot power for both hemispheres and their sum, and highlight the optimal time point
        subplot(ceil(numSubjects / 6), 6, subj);
%         plot(timePoints, powerLeft(subj, :), 'b', 'LineWidth', 1); % Plot left hemisphere power as blue
        hold on;
%         plot(timePoints, powerRight(subj, :), 'r', 'LineWidth', 1); % Plot right hemisphere power as red
        plot(timePoints, combinedPower, 'k', 'LineWidth', 1.2); % Plot combined power as magenta
        if ~isnan(optimalIdx)
            plot(timePoints(optimalIdx), optimalPower, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k'); % Highlight the optimal point
        end
        hold off;
        axis tight;
        title(sprintf('Subject %s', ID{subj}));
        title(sprintf('Subject %s', num2str((subj))));
    end
end

if plot_flag
    % Additional plotting settings if required
    legend({'Left Hemisphere Power', 'Right Hemisphere Power', 'Optimal Time Point'}, 'Location', 'best');
    set(gcf, 'Position', [100, 100, 1200, 1000]);
end

% Calculate the histogram of optimal time points if needed to find the most common interval
if nargout > 1
    figure;
    h = histogram(optimalTimePoints, 'FaceColor', 'k');
    title('Distribution of Optimal Time Points Across Subjects');
    xlabel('Time Points');
    ylabel('Frequency');
    
    % Calculate an interval covering 75% of the peak data
    binEdges = h.BinEdges;
    binCounts = h.Values;
    [maxCount, idx] = max(binCounts);
    cumCounts = cumsum(binCounts);
    total = cumCounts(end);
    threshold = 0.75 * total;
    lowerIdx = find(cumCounts > threshold, 1, 'first');
    interval = [binEdges(idx), binEdges(lowerIdx + 1)];
end



end
