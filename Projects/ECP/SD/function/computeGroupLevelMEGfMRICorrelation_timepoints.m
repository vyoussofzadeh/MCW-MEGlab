function [groupCorrelation, optimalTimePoints] = computeGroupLevelMEGfMRICorrelation_timepoints(MEG_LI, fMRI_LI, timePoints, lowerBound, upperBound)
    % Find individual optimal time points
    optimalTimePoints = findIndividualOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, NaN, lowerBound, upperBound); % NaN means no plot
    
    numSubjects = size(MEG_LI, 1);
    optimalMEG_LI = zeros(numSubjects, 1);
    
    % Extract MEG LI at optimal time points for each subject
    for subj = 1:numSubjects
        timePointIndex = timePoints == optimalTimePoints(subj);
        optimalMEG_LI(subj) = MEG_LI(subj, timePointIndex);
    end

    % Perform group-level correlation analysis
    validIndices = ~isnan(optimalMEG_LI) & ~isnan(fMRI_LI); % Exclude NaN values
    groupCorrelation = corr(optimalMEG_LI(validIndices), fMRI_LI(validIndices), 'Rows', 'complete');
    
    % Display or return the correlation result
%     disp(['Group-level correlation between MEG LI at optimal time points and fMRI LI: ', num2str(groupCorrelation)]);
end
