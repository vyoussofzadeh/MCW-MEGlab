function concordance = calculateConcordance(MEG_LI, fMRI_LI, timePoints, optimalIntervals)
    numSubjects = size(MEG_LI, 1);
    concordantCount = 0;

    for subj = 1:numSubjects
        % Calculate mean MEG LI within the optimal interval
        intervalStartIdx = find(timePoints == optimalIntervals(subj, 1));
        intervalEndIdx = find(timePoints == optimalIntervals(subj, 2));
        meanOptimalMEG_LI = mean(MEG_LI(subj, intervalStartIdx:intervalEndIdx));

        % Classify MEG LI
        if meanOptimalMEG_LI > 0.20
            MEG_Class = 1;
        elseif meanOptimalMEG_LI < -0.20
            MEG_Class = -1;
        else
            MEG_Class = 0;
        end

        % Classify fMRI LI
        if fMRI_LI(subj) > 0.10
            fMRI_Class = 1;
        elseif fMRI_LI(subj) < -0.10
            fMRI_Class = -1;
        else
            fMRI_Class = 0;
        end

        % Check Concordance
        if MEG_Class == fMRI_Class
            concordantCount = concordantCount + 1;
        end
    end

    % Calculate Concordance Percentage
    concordance = (concordantCount / numSubjects) * 100;
end
