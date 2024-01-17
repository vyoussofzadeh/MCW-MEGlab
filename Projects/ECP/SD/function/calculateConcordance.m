function [concordance, discordantSubs] = calculateConcordance(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre,timePoints, optimalIntervals)
numSubjects = size(MEG_LI, 1);
concordantCount = 0;
discordantSubs = []; % Initialize an empty array to store indices of discordant subjects

for subj = 1:numSubjects
    % Calculate mean MEG LI within the optimal interval
    intervalStartIdx = find(timePoints == optimalIntervals(subj, 1));
    intervalEndIdx = find(timePoints == optimalIntervals(subj, 2));
    meanOptimalMEG_LI = mean(MEG_LI(subj, intervalStartIdx:intervalEndIdx));
    
    % Classify MEG LI
    if meanOptimalMEG_LI > MEG_thre
        MEG_Class = 1;
    elseif meanOptimalMEG_LI < -1*MEG_thre
        MEG_Class = -1;
    else
        MEG_Class = 0;
    end
    
    % Classify fMRI LI
    if fMRI_LI(subj) > fMRI_thre
        fMRI_Class = 1;
    elseif fMRI_LI(subj) < -1*fMRI_thre
        fMRI_Class = -1;
    else
        fMRI_Class = 0;
    end
    
    % Check Concordance
    if MEG_Class == fMRI_Class
        concordantCount = concordantCount + 1;
    else
        % Store the index of discordant subjects
        discordantSubs = [discordantSubs, subj];
    end
end

% Calculate Concordance Percentage
concordance = (concordantCount / numSubjects) * 100;
end
