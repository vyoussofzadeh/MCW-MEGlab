function [optimalInterval, correlations] = findOptimalMEGInterval(MEG_LI, fMRI_LI, timePoints, intervalSize)
    
    % Initialize variables
    numIntervals = length(timePoints) - intervalSize + 1;
    correlations = zeros(1, numIntervals);
    
    % Loop through each interval
    for i = 1:numIntervals
        intervalMEG_LI = mean(MEG_LI(:, i:(i + intervalSize - 1)), 2);
        
        % Compute correlation for each interval
        corrMatrix = corrcoef(intervalMEG_LI, fMRI_LI);
        correlations(i) = corrMatrix(1, 2);
    end
    
    % Identify the interval with the highest correlation
    [~, maxIdx] = max(correlations);
    optimalInterval = timePoints(maxIdx:maxIdx + intervalSize - 1);
    
    % Optional: Plot the correlation values across intervals
    figure;
    plot(timePoints, correlations);
    xlabel('Interval Start Time (sec)');
    ylabel('Correlation with fMRI LI');
    title('Correlation of MEG LI with fMRI LI across Intervals');
end

% Usage
% optimalInterval, correlations = findOptimalMEGInterval(MEG_LI_Data, fMRI_LI_Data, TimePoints, IntervalSize);
