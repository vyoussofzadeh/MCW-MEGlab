function [groupCorrelation, optimalIntervals] = computeGroupLevelMEGfMRICorrelation(MEG_LI, fMRI_LI, timePoints, intervalSize, stepSize)
    % Find individual optimal intervals
    optimalIntervals = findIndividualOptimalIntervals(MEG_LI, fMRI_LI, timePoints, intervalSize, stepSize, NaN); % NaN means no plot
    
    numSubjects = size(MEG_LI, 1);
    meanOptimalMEG_LI = zeros(numSubjects, 1);
    
    % Compute mean MEG LI within optimal intervals for each subject
    for subj = 1:numSubjects
        intervalStart = optimalIntervals(subj, 1);
        intervalEnd = optimalIntervals(subj, 2);
        intervalIndex = timePoints >= intervalStart & timePoints <= intervalEnd;
        meanOptimalMEG_LI(subj) = mean(MEG_LI(subj, intervalIndex));
    end

    % Perform group-level correlation analysis
    validIndices = ~isnan(meanOptimalMEG_LI) & ~isnan(fMRI_LI); % Exclude NaN values
    groupCorrelation = corr(meanOptimalMEG_LI(validIndices), fMRI_LI(validIndices), 'Rows', 'complete');
    
    % Display or return the correlation result
    disp(['Group-level correlation between mean MEG LI and fMRI LI: ', num2str(groupCorrelation)]);
end

% Usage Example:
% [groupCorrelation, optimalIntervals] = computeGroupLevelMEGfMRICorrelation(MEG_LI_Data, fMRI_LI_Data, TimePoints, IntervalSize, StepSize);
