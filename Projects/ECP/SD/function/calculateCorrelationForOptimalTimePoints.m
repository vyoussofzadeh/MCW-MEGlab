function [correlationCoeff, pValue, optimalMEGLI, optimalfMRILI] = calculateCorrelationForOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, optimalTimePoints)
% Calculates the correlation coefficient between MEG and fMRI LI values
% at optimal time points identified for each subject.

numSubjects = size(MEG_LI, 1);
optimalMEGLI = zeros(numSubjects, 1);
optimalfMRILI = zeros(numSubjects, 1);

for subj = 1:numSubjects
    % Extract MEG LI at the optimal time point for each subject
    optimalTimePointIdx = find(timePoints == optimalTimePoints(subj), 1);
    if ~isempty(optimalTimePointIdx)
        optimalMEGLI(subj) = MEG_LI(subj, optimalTimePointIdx);
    else
        optimalMEGLI(subj) = NaN; % Handle cases where the optimal time point is not found
    end
    
    % Extract the corresponding fMRI LI for each subject
    optimalfMRILI(subj) = fMRI_LI(subj);
end

% Remove subjects with NaN values from the analysis (if any)
validIdx = ~isnan(optimalMEGLI) & ~isnan(optimalfMRILI);
optimalMEGLI = optimalMEGLI(validIdx);
optimalfMRILI = optimalfMRILI(validIdx);

% Calculate Pearson correlation coefficient and p-value
[correlationCoeff, pValue] = corr(optimalMEGLI, optimalfMRILI);

figure, bar([optimalMEGLI,optimalfMRILI]), 
lgnd = legend({'MEG', 'fMRI'});
set(gcf, 'Position', [100   200   1500   300]);
set(gca,'color','none');
set(lgnd,'color','none');

disp(['Pearson correlation coefficient between optimal MEG LI and fMRI LI: ', num2str(correlationCoeff)]);
disp(['P-value: ', num2str(pValue)]);

end
