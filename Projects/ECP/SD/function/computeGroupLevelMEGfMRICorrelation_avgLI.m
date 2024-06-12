function [groupCorrelation] = computeGroupLevelMEGfMRICorrelation_avgLI(MEG_LI, fMRI_LI, timePoints, lowerBound, upperBound)
    % Validate input dimensions
    assert(size(MEG_LI,2) == length(timePoints), 'Dimension mismatch between MEG_LI and timePoints');
    
    % Find indices of time points within specified bounds
    validTimeIndices = timePoints >= lowerBound & timePoints <= upperBound;
    
    numSubjects = size(MEG_LI, 1);
    avgMEG_LI = zeros(numSubjects, 1);
    
    % Compute average MEG LI within bounds for each subject
    for subj = 1:numSubjects
        validMEG_LIs = MEG_LI(subj, validTimeIndices);
%         avgMEG_LI(subj) = mean(validMEG_LIs, 'omitnan');
        avgMEG_LI(subj) = sum(validMEG_LIs, 'omitnan');
    end

    % Perform group-level correlation analysis
    validIndices = ~isnan(avgMEG_LI) & ~isnan(fMRI_LI); % Exclude NaN values
    groupCorrelation = corr(avgMEG_LI(validIndices), fMRI_LI(validIndices), 'Rows', 'complete');
    
    % Display or return the correlation result
    disp(['Group-level correlation between average MEG LI and fMRI LI: ', num2str(groupCorrelation)]);
end

