function [groupCorrelation, optimalInterval, optimalInterval_all, pval] = computeGroupLevelMEGfMRICorrelation_timepoints_interval_rSNR(MEG_LI,rSNR_left, rSNR_right, fMRI_LI, timePoints, lowerBound, upperBound)
    % Find individual optimal time points
%     optimalTimePoints = findIndividualOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, NaN, lowerBound, upperBound); % NaN means no plot
    [optimalIndices, optimalMEG_LI] = findIndividualOptimalTimePoints_interval_rSNR(rSNR_left, rSNR_right, fMRI_LI, timePoints, NaN, lowerBound, upperBound); % NaN means no plot

    %     numSubjects = size(MEG_LI, 1);
%     optimalMEG_LI = zeros(numSubjects, 1);
    
    optimalInterval = timePoints(optimalIndices,:);
    
%     % Extract MEG LI at optimal time points for each subject
%     for subj = 1:numSubjects
%         optimalInterval(subj,:) = timePoints(optimalIndices(subj),:);
%         timePointIndex = timePoints == optimalTimePoints(subj);
%         optimalMEG_LI(subj) = MEG_LI(subj, timePointIndex);
%     end

    % Perform group-level correlation analysis
%     validIndices = ~isnan(optimalMEG_LI) & ~isnan(fMRI_LI); % Exclude NaN values
    [groupCorrelation, pval] = corr(optimalMEG_LI, fMRI_LI, 'Rows', 'complete');
    
    % Display or return the correlation result
    disp(['Group-level correlation between MEG LI at optimal time points and fMRI LI: ', num2str(groupCorrelation)]);
    
    
    %%
    optimalInterval_all = timePoints(mode(optimalIndices),:);
    
end
