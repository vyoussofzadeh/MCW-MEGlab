function [optimalTimePoints, interval] = plotSubjectPowerOverlay(powerLeft, powerRight, ID, timePoints, plot_flag, lowerBound, upperBound)
% Plots overlaid power values for left and right hemispheres for all subjects with individual peak annotations.
% Searches for the optimal time point based on the maximum absolute power difference strictly within the interval defined by lowerBound and upperBound.

numSubjects = size(powerLeft, 1);
optimalTimePoints = zeros(numSubjects, 1); % Initialize array to store optimal time points

% Create a figure to hold all subplots if plotting is enabled
if plot_flag == 1, figure; end

for subj = 1: numSubjects
    maxPowerDiff = 0; % Initialize maximum absolute MEG LI
    optimalIdx = NaN; % Initialize index of optimal time point
    
    %     % Calculate the difference in power between left and right hemispheres for plotting
    powerDiff = powerLeft(subj, :) - powerRight(subj, :);
    
    for i = 1:length(timePoints)
        % Check if the current time point is within the desired interval
        if timePoints(i) >= lowerBound && timePoints(i) <= upperBound
            % Find the maximum absolute MEG LI within the interval
            if abs(powerDiff(i)) > maxPowerDiff
                maxPowerDiff = abs(powerDiff(i));
                optimalIdx = i;
            end
        end
    end
    
    
    if ~isnan(optimalIdx)
        optimalTimePoints(subj) = timePoints(optimalIdx);
    else
        optimalTimePoints(subj) = NaN; % Assign NaN if no maximum is found in the interval
    end
    
    if plot_flag == 1
        subplot(ceil(numSubjects / 6), 6, subj); % Adjust layout to fit all subjects
        plot(timePoints, powerLeft(subj, :), 'b'); % Plot left hemisphere power as a blue line
        hold on;
        plot(timePoints, powerRight(subj, :), 'r'); % Plot right hemisphere power as a red line
        plot(timePoints, powerDiff, 'k', 'LineWidth', 1.5); % Overlay the difference as a dashed black line
        % Highlight the optimal time point, if within bounds
        if optimalIdx > 0
            plot(optimalTimePoints(subj), powerDiff(optimalIdx), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'y'); % Highlight the optimal time point
        end
        hold off;
        axis tight;
        title(sprintf('Subject %s', ID{subj}));
        title(sprintf('Subject %s', num2str((subj))));
    end
    %     ylim([0,3])
end

interval = [];

if plot_flag == 1
    xlabel('Time Points');
    ylabel('Power (dB)');
    lgd = legend({'Left Power', 'Right Power', 'Left-Right Diff', 'Peak'}, 'Location', 'best', 'Orientation', 'horizontal');
    lgdPos = lgd.Position; % Get current position
    lgdPos(2) = lgdPos(2) - 0.11; % Move legend down
    lgd.Position = lgdPos; % Set new position
    set(gcf, 'Position', [100   100   1200   1000]);
    
    xlabel('Time Points');
    ylabel('Power (dB)');
    lgd = legend({'Left Power', 'Right Power', 'Optimal Time Point'}, 'Location', 'best', 'Orientation', 'horizontal');
    lgdPos = lgd.Position; % Get current position
    lgdPos(2) = lgdPos(2) - 0.11; % Move legend down
    lgd.Position = lgdPos; % Set new position
    set(gcf, 'Position', [100, 100, 1200, 1000]);
    
%     % Plot distribution of peaks for left and right power source values
%     figure;
%     subplot(1, 2, 1);
%     histogram(peakLeft, 'FaceColor', 'b');
%     title({'Distribution of'; 'Left Hemisphere Peak Source Mag'});
%     xlabel('Peak Power (dB)');
%     ylabel('Number of Subjects');
%     
%     subplot(1, 2, 2);
%     histogram(peakRight, 'FaceColor', 'r');
%     title({'Distribution of'; 'Right Hemisphere Peak Source Mag'});
%     xlabel('Peak Power (dB)');
%     ylabel('Number of Subjects');
    
    % Plot distribution of latency peaks across subjects
    figure;
    %     histogram(optimalTimePoints, 'FaceColor', 'k');
    h = histogram(optimalTimePoints, 'FaceColor', 'k');
    text(mean(optimalTimePoints), max(h.Values), sprintf('Mean: %.2f', mean(optimalTimePoints)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');
    
    title('Distribution of Latency Peaks');
    xlabel('Time Points');
    ylabel('Number of Subjects');
    
    %     h = histogram(optimalTimePoints, 'FaceColor', 'k');
    binEdges = h.BinEdges;
    binCounts = h.Values;
    
    % Example: Find the range that contains 75% of the data starting from the mode
    [maxCount, idx] = max(binCounts); % Find the bin with the most data
    cumCounts = cumsum(binCounts);
    total = cumCounts(end);
    
    % Define the cumulative percentage threshold
    threshold = 0.75 * total;
    lowerBound = binEdges(idx);
    upperBound = binEdges(find(cumCounts > threshold, 1, 'first') + 1);
    
    fprintf('Approximately 75%% of the data is between %.2f and %.2f.\n', lowerBound, upperBound);
    
    
    interval = [lowerBound, upperBound];
     
end
end


