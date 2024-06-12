function plotOptimalTimePointsOnMEG(MEG_LI, fMRI_LI, ID, timePoints, optimalTimePoints)

numSubjects = size(MEG_LI, 1);

% Create a figure to hold all subplots
figure;

for subj = 1:numSubjects
    % Find the index of the optimal time point
    optimalTimePointIdx = find(timePoints == optimalTimePoints(subj));
    
    subplot(ceil(numSubjects / 6), 6, subj); % Adjust for 5 columns
    plot(timePoints, MEG_LI(subj, :), 'b'); % Plot MEG LI as a blue line
    hold on;
    
    % Highlight the optimal time point
    plot(timePoints(optimalTimePointIdx), MEG_LI(subj, optimalTimePointIdx), 'ro'); % Mark optimal point as a red circle
    
    % Set the y-axis limits
    ylim([-100 100]);
    
    hold off;
%     axis tight
    set(gca,'color','none');
%     title(sprintf('%s | fMRI LI: %.2f | MEG LI: %.2f', ID{subj}, fMRI_LI(subj), MEG_LI(subj, optimalTimePointIdx)));
    title(sprintf('S %d | fMRI LI: %.2f | MEG LI: %.2f', subj, fMRI_LI(subj), MEG_LI(subj, optimalTimePointIdx)));

end
xlabel('Time Points');
ylabel('MEG LI');
set(gcf, 'Position', [100   100   1200   1000]);
% set(gcf, 'Position', [100   100   1800   1500]);
end
