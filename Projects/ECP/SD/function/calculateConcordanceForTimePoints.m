function [concordance, discordantSubs, optimalMEG_LI] = calculateConcordanceForTimePoints(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimePoints)

numSubjects = size(MEG_LI, 1);
concordantCount = 0;
discordantSubs = []; % Initialize an empty array to store indices of discordant subjects

for subj = 1:numSubjects
    
    % Extract MEG LI at the optimal time point
    optimalTimePointIdx = timePoints == optimalTimePoints(subj);
    optimalMEG_LI(subj) = MEG_LI(subj, optimalTimePointIdx);
%     optimalMEG_LI = optimalMEG_LI./max(optimalMEG_LI(:));
    
    % Classify MEG LI
    if optimalMEG_LI(subj) > MEG_thre
        MEG_Class = 1;
    elseif optimalMEG_LI < -MEG_thre
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
%         LI_concordance = [LI_concordance,optimalMEG_LI(subj)];
    end
end

% Calculate Concordance Percentage
concordance = (concordantCount / numSubjects) * 100;

disp(['Group-level concordance between MEG LI at optimal time points and fMRI LI: ', num2str(concordance)]);

end



