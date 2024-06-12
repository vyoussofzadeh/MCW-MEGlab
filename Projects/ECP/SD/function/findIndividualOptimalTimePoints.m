
function [optimalTimePoints, LI_opt] = findIndividualOptimalTimePoints(MEG_LI, fMRI_LI, timePoints, subjectsForPlot, lowerBound, upperBound)
numSubjects = size(MEG_LI, 1);
optimalTimePoints = zeros(numSubjects, 1); % Initialize array for optimal time points

% Create a figure outside the loop for all subplots
if ~isempty(subjectsForPlot) && mean(~isnan(subjectsForPlot))
    figure;
end

k=1;
for subj = 1:numSubjects
    maxMEG_LI = 0; % Initialize maximum absolute MEG LI
    optimalTimePointIdx = NaN; % Initialize index of optimal time point
    
    for i = 1:length(timePoints)
        % Check if the current time point is within the desired interval
        if timePoints(i) >= lowerBound && timePoints(i) <= upperBound
            % Find the maximum absolute MEG LI within the interval
            if abs(MEG_LI(subj, i)) > maxMEG_LI
                maxMEG_LI = abs(MEG_LI(subj, i));
                optimalTimePointIdx = i;
            end
        end
    end
    
    if ~isnan(optimalTimePointIdx)
        optimalTimePoints(subj) = timePoints(optimalTimePointIdx);
    else
        optimalTimePoints(subj) = NaN; % Assign NaN if no maximum is found in the interval
    end
    
    % Optional Plotting for selected subjects
    if ismember(subj, subjectsForPlot)
        % Create subplot for each selected subject
        subplot(ceil(length(subjectsForPlot) / 2), 2, find(subjectsForPlot == subj)); % Adjust for 2 columns layout
        plot(timePoints, MEG_LI(subj, :), 'b'); % Plot MEG LI as a blue line
        hold on;
        
        % Highlight the optimal time point
        plot(optimalTimePoints(subj), MEG_LI(subj, optimalTimePointIdx), 'ro'); % Mark optimal point as a red circle
        
        % Set the y-axis limits
        ylim([-100 100]);
        
        hold off;
        title(sprintf('Subject %d | fMRI LI: %.2f | MEG LI: %.2f', subj, fMRI_LI(subj), MEG_LI(subj, optimalTimePointIdx)));
        xlabel('Time Points');
        ylabel('MEG LI');
        set(gcf, 'Position', [200   400   800   400]);

        
        LI_opt(k) = MEG_LI(subj, optimalTimePointIdx);
        k = k+1;
    end
end
end

