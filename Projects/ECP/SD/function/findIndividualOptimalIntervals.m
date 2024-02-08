function optimalIntervals = findIndividualOptimalIntervals(MEG_LI, fMRI_LI, timePoints, intervalSize, stepSize, subjectForPlot)
numSubjects = size(MEG_LI, 1);
numIntervals = floor((length(timePoints) - intervalSize) / stepSize) + 1;
optimalIntervals = zeros(numSubjects, 2); % Start and end of optimal interval

for subj = 1:numSubjects
    minDiff = inf; % Initialize minimum difference
    firstMinDiffIdx = NaN; % Initialize index of first minimum difference
    
    for i = 1:numIntervals
        intervalStartIdx = ((i - 1) * stepSize) + 1;
        intervalEndIdx = intervalStartIdx + intervalSize - 1;
        
        % Ensure the interval does not exceed the bounds of timePoints
        if intervalEndIdx > length(timePoints)
            intervalEndIdx = length(timePoints);
            intervalStartIdx = intervalEndIdx - intervalSize + 1;
        end
        
        intervalMEG_LI = MEG_LI(subj, intervalStartIdx:intervalEndIdx);
        meanIntervalMEG_LI = mean(intervalMEG_LI); % Mean MEG LI for the interval
        
        % Check the difference with fMRI LI
        diff = abs(meanIntervalMEG_LI - fMRI_LI(subj));
        if diff < minDiff
            minDiff = diff;
            firstMinDiffIdx = i;
        end
    end
    
    optimalStartIdx = ((firstMinDiffIdx - 1) * stepSize) + 1;
    optimalEndIdx = optimalStartIdx + intervalSize - 1;
    
    % Ensure the optimal interval does not exceed the bounds of timePoints
    if optimalEndIdx > length(timePoints)
        optimalEndIdx = length(timePoints);
    end
    
    optimalIntervals(subj, :) = [timePoints(optimalStartIdx), timePoints(optimalEndIdx)];
    
    % Optional Plotting
    if subj == subjectForPlot
        figure,
        intervalStartIdx = find(timePoints == optimalIntervals(subj, 1));
        intervalEndIdx = find(timePoints == optimalIntervals(subj, 2));
        meanOptimalMEG_LI = mean(MEG_LI(subj, intervalStartIdx:intervalEndIdx));
        
        plot(timePoints, MEG_LI(subj, :), 'b'); % Plot MEG LI as a blue line
        hold on;
        
        % Highlight the optimal interval
        patch([timePoints(intervalStartIdx), timePoints(intervalEndIdx), timePoints(intervalEndIdx), timePoints(intervalStartIdx)], ...
            [min(MEG_LI(subj, :)), min(MEG_LI(subj, :)), max(MEG_LI(subj, :)), max(MEG_LI(subj, :))], ...
            'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Adjust color and transparency as needed
        
        % Set the y-axis limits
        ylim([-100 100]);
        
        hold off;
        title(sprintf('Subject %d | fMRI LI: %.2f | Mean MEG LI: %.2f', subj, fMRI_LI(subj), meanOptimalMEG_LI));
        xlabel('Time Points');
        ylabel('MEG LI');
    end
end