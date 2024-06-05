function [concordance, discordantSubs, groupCorrelation] = calculateConcordanceForTimePoints_interval(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimeInterval)

numSubjects = size(MEG_LI, 1);
concordantCount = 0;
discordantSubs = []; % Initialize an empty array to store indices of discordant subjects

% MEG_thre = MEG_thre*100;
% fMRI_thre = fMRI_thre*100;

for subj = 1:numSubjects
    
    % Find the indices within the optimal time interval
    intervalIndices = abs(timePoints(:,1) - optimalTimeInterval(subj,1)) < 0.01 & abs(timePoints(:,2) - optimalTimeInterval(subj,2)) < 0.01;
    
    % Extract MEG LI within the optimal time interval and find the maximum absolute value
%     [~, maxIdx] = max(abs(MEG_LI(subj, timePoints)));
%     optimalMEG_LI = MEG_LI(subj, intervalIndices(maxIdx));
    
    optimalMEG_LI(subj) = MEG_LI(subj,intervalIndices);
%     optimalMEG_LI = optimalMEG_LI./max(optimalMEG_LI(:));
    
    % Classify MEG LI
    if optimalMEG_LI(subj) > MEG_thre
        MEG_Class = 1;
    elseif optimalMEG_LI(subj) < -MEG_thre
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

groupCorrelation = corr(optimalMEG_LI', fMRI_LI, 'Rows', 'complete');

disp(['Group-level concordance between MEG LI at optimal time points and fMRI LI: ', num2str(concordance), '%']);

end
