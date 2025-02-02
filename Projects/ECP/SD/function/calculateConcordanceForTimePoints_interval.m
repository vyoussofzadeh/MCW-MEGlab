function [concordance, discordantSubs, groupCorrelation, optimalMEG_LI, Corrpval, kappa] = calculateConcordanceForTimePoints_interval(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimeInterval)

numSubjects = size(MEG_LI, 1);
concordantCount = 0;
discordantSubs = []; % Initialize an empty array to store indices of discordant subjects
optimalMEG_LI = zeros(numSubjects, 1);

% Preallocate classification vectors
MEG_Class = zeros(numSubjects, 1);
fMRI_Class = zeros(numSubjects, 1);

% 2) If optimalTimeInterval is a single row (1 x 2), replicate it for all subjects
if size(optimalTimeInterval,1) == 1 && size(optimalTimeInterval,2) == 2
    optimalTimeInterval = repmat(optimalTimeInterval, numSubjects, 1);
elseif size(optimalTimeInterval,1) ~= numSubjects
    error(['optimalTimeInterval must either be (#Subjects x 2) or (1 x 2). ' ...
        'Currently it is (%d x %d).'], size(optimalTimeInterval,1), size(optimalTimeInterval,2));
end

for subj = 1:numSubjects
    % Assuming optimalTimeInterval contains exact start and stop times for each subject
    intervalIndices = timePoints(:,1) >= optimalTimeInterval(subj, 1) & timePoints(:,2) <= optimalTimeInterval(subj, 2);
    if any(intervalIndices)
        % Extract the range of MEG LI values within the optimal time interval
        subjectMEG_LIs = MEG_LI(subj, intervalIndices);
        % Take the average or maximum based on your requirement
        optimalMEG_LI(subj) = mean(subjectMEG_LIs);  % Example use of mean
        
        % Classify MEG LI
        if optimalMEG_LI(subj) > MEG_thre
            MEG_Class(subj) = 1;
        elseif optimalMEG_LI(subj) < -MEG_thre
            MEG_Class(subj) = -1;
        else
            MEG_Class(subj) = 0;
        end
        
        % Classify fMRI LI
        if fMRI_LI(subj) > fMRI_thre
            fMRI_Class(subj) = 1;
        elseif fMRI_LI(subj) < -fMRI_thre
            fMRI_Class(subj) = -1;
        else
            fMRI_Class(subj) = 0;
        end
        
        % Check Concordance
        if MEG_Class(subj) == fMRI_Class(subj)
            concordantCount = concordantCount + 1;
        else
            % Store the index of discordant subjects
            discordantSubs = [discordantSubs, subj];
        end
    else
        % Handle cases where no time points are found within the interval
        optimalMEG_LI(subj) = NaN;  % Optional: handle as you see fit
    end
end

% Calculate Concordance Percentage
% concordance = (concordantCount / numSubjects) * 100;
% 
% % Calculate correlation between MEG LI and fMRI LI across subjects
% % [groupCorrelation, pval] = corr(optimalMEG_LI, fMRI_LI, 'Rows', 'complete');
% [groupCorrelation, pval] = corr(optimalMEG_LI, fMRI_LI, 'Type', 'Spearman');
% 
% cfg = []; cfg.thre = fMRI_thre; cfg.LI = fMRI_LI; fmri_LIs_trn = do_ternary_classification2(cfg);
% cfg = []; cfg.thre = MEG_thre; cfg.LI = optimalMEG_LI; MEG_LIs_trn = do_ternary_classification2(cfg);
% 
% 
% confMat = confusionmat(MEG_LIs_trn, fmri_LIs_trn);
% 
% % Suppose confMat is your confusion matrix (square matrix of size NxN).
% totalObservations = sum(confMat(:));
% correctMatches = trace(confMat);  % sum of diagonal
% concordance = correctMatches / totalObservations;
% disp(['Concordance (overall agreement) = ' num2str(concordance)]);

% Calculate correlation between MEG LI and fMRI LI across subjects
[groupCorrelation, Corrpval] = corr(optimalMEG_LI, fMRI_LI, 'Type', 'Spearman');

cfg = []; cfg.thre = fMRI_thre; cfg.LI = fMRI_LI; fmri_LIs_trn = do_ternary_classification2(cfg);
cfg = []; cfg.thre = MEG_thre; cfg.LI = optimalMEG_LI; MEG_LIs_trn = do_ternary_classification2(cfg);

%%
confMat = confusionmat(MEG_LIs_trn, fmri_LIs_trn);

% Suppose confMat is your confusion matrix (square matrix of size NxN).
totalObservations = sum(confMat(:));
correctMatches = trace(confMat);  % sum of diagonal
concordance = correctMatches / totalObservations * 100;
disp(['Concordance (overall agreement) = ' num2str(concordance)]);

[kappa, po, pe] = ecpfunc_cohenKappa(confMat);
disp(['Cohen''s kappa = ' num2str(kappa)]);



end

% function [concordance, discordantSubs, groupCorrelation, pval] = calculateConcordanceForTimePoints_interval(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimeInterval)
%
% numSubjects = size(MEG_LI, 1);
% concordantCount = 0;
% discordantSubs = []; % Initialize an empty array to store indices of discordant subjects
% optimalMEG_LI = zeros(numSubjects, 1);
% MEG_Class = zeros(numSubjects, 1);
% fMRI_Class = zeros(numSubjects, 1);
%
% for subj = 1:numSubjects
%     % Assuming optimalTimeInterval contains exact start and stop times for each subject
%     intervalIndices = timePoints(:,1) >= optimalTimeInterval(subj, 1) & timePoints(:,2) <= optimalTimeInterval(subj, 2);
%     if any(intervalIndices)
%         % Extract the range of MEG LI values within the optimal time interval
%         subjectMEG_LIs = MEG_LI(subj, intervalIndices);
%         % Aggregate the MEG LI values within the interval, e.g., by taking the mean
%         optimalMEG_LI(subj) = mean(subjectMEG_LIs);  % Example use of mean
%
%         % Classify MEG LI
%         if optimalMEG_LI(subj) > MEG_thre
%             MEG_Class(subj) = 1;
%         elseif optimalMEG_LI(subj) < -MEG_thre
%             MEG_Class(subj) = -1;
%         else
%             MEG_Class(subj) = 0;
%         end
%
%         % Classify fMRI LI
%         if fMRI_LI(subj) > fMRI_thre
%             fMRI_Class(subj) = 1;
%         elseif fMRI_LI(subj) < -fMRI_thre
%             fMRI_Class(subj) = -1;
%         else
%             fMRI_Class(subj) = 0;
%         end
%
%         % Check Concordance
%         if MEG_Class(subj) == fMRI_Class(subj)
%             concordantCount = concordantCount + 1;
%         else
%             discordantSubs = [discordantSubs, subj];  % Record discordant subject indices
%         end
%     end
% end
%
% % Calculate Concordance Percentage
% concordance = (concordantCount / numSubjects) * 100;
%
% % Calculate correlation between classified MEG LI and fMRI LI across subjects
% [groupCorrelation, pval] = corr(MEG_Class, fMRI_Class, 'Rows', 'complete', 'Type', 'Spearman');  % Using Spearman as it is more appropriate for ordinal data
%
% end


% function [concordance, discordantSubs, groupCorrelation, pval] = calculateConcordanceForTimePoints_interval(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, timePoints, optimalTimeInterval)
%
% numSubjects = size(MEG_LI, 1);
% concordantCount = 0;
% discordantSubs = []; % Initialize an empty array to store indices of discordant subjects
%
% % MEG_thre = MEG_thre*100;
% % fMRI_thre = fMRI_thre*100;
%
% for subj = 1:numSubjects
%
%     % Find the indices within the optimal time interval
%     intervalIndices = abs(timePoints(:,1) - optimalTimeInterval(subj,1)) < 0.001 & abs(timePoints(:,2) - optimalTimeInterval(subj,2)) < 0.001;
%
%     % Extract MEG LI within the optimal time interval and find the maximum absolute value
% %     [~, maxIdx] = max(abs(MEG_LI(subj, timePoints)));
% %     optimalMEG_LI = MEG_LI(subj, intervalIndices(maxIdx));
%
%     optimalMEG_LI(subj) = MEG_LI(subj,intervalIndices);
% %     optimalMEG_LI = optimalMEG_LI./max(optimalMEG_LI(:));
%
%     % Classify MEG LI
%     if optimalMEG_LI(subj) > MEG_thre
%         MEG_Class = 1;
%     elseif optimalMEG_LI(subj) < -MEG_thre
%         MEG_Class = -1;
%     else
%         MEG_Class = 0;
%     end
%
%     % Classify fMRI LI
%     if fMRI_LI(subj) > fMRI_thre
%         fMRI_Class = 1;
%     elseif fMRI_LI(subj) < -fMRI_thre
%         fMRI_Class = -1;
%     else
%         fMRI_Class = 0;
%     end
%
%     % Check Concordance
%     if MEG_Class == fMRI_Class
%         concordantCount = concordantCount + 1;
%     else
%         % Store the index of discordant subjects
%         discordantSubs = [discordantSubs, subj];
%     end
% end
%
% % Calculate Concordance Percentage
% concordance = (concordantCount / numSubjects) * 100;
%
% [groupCorrelation, pval] = corr(optimalMEG_LI', fMRI_LI, 'Rows', 'complete');
% % groupCorrelation = corr2(optimalMEG_LI', fMRI_LI);
%
%
% % disp(['Group-level concordance between MEG LI at optimal time points and fMRI LI: ', num2str(concordance), '%']);
%
% end
