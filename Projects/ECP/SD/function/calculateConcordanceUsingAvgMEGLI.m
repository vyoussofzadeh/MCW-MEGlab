function [concordance, discordantSubs] = calculateConcordanceUsingAvgMEGLI(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, lowerBound, upperBound)

numSubjects = size(MEG_LI, 1);
concordantCount = 0;
discordantSubs = []; % Initialize an empty array to store indices of discordant subjects

% Find indices of time points within specified bounds
validTimeIndices = timePoints >= lowerBound & timePoints <= upperBound;

for subj = 1:numSubjects
    
    % Compute average MEG LI within bounds
    validMEG_LIs = MEG_LI(subj, validTimeIndices);
    avgMEG_LI = mean(validMEG_LIs, 'omitnan');
    
    % Classify average MEG LI
    if avgMEG_LI > MEG_thre
        MEG_Class = 1;
    elseif avgMEG_LI < -MEG_thre
        MEG_Class = -1;
    else
        MEG_Class = 0;
    end
    
    % Classify fMRI LI
    if fMRI_LI(subj) > fMRI_thre
        fMRI_Class = 1;
    elseif fMRI_LI(subj) < -fMRI_thre
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

disp(['Group-level concordance between average MEG LI within bounds and fMRI LI: ', num2str(concordance)]);

end
