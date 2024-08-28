function [optimalIndices, maxLI_opt] = findIndividualOptimalTimePoints_interval_peak(MEG_LI, fMRI_LI, timePoints, subjectsForPlot, lowerBound, upperBound)
numSubjects = size(MEG_LI, 1);
optimalIndices = zeros(numSubjects, 1); % Initialize array for optimal indices
maxLI_opt = zeros(numSubjects, 1); % Initialize array for maximum LI

% Create a figure outside the loop for all subplots
if ~isempty(subjectsForPlot) && mean(~isnan(subjectsForPlot))
    figure;
end

for subj = 1:numSubjects
    localPeaks = []; % Initialize array for local peaks
    peakIndices = []; % Initialize array for indices of local peaks
    
    for i = 2:length(timePoints)-1
        % Check if the current time point is within the desired interval
        if timePoints(i, 1) >= lowerBound && timePoints(i, 2) <= upperBound
            % Find local peaks
            if MEG_LI(subj, i) > MEG_LI(subj, i-1) && MEG_LI(subj, i) > MEG_LI(subj, i+1)
                localPeaks = [localPeaks, MEG_LI(subj, i)];
                peakIndices = [peakIndices, i];
            end
        end        
    end
    
    if ~isempty(localPeaks)
        [maxLI, maxIdx] = max(localPeaks);
        optimalTimePointIdx = peakIndices(maxIdx);
        maxLI_opt(subj) = maxLI;
    else
        optimalTimePointIdx = NaN;
        maxLI_opt(subj) = NaN;
    end
    
    if ~isnan(optimalTimePointIdx)
        optimalIndices(subj) = optimalTimePointIdx; % Store the index of the optimal time point
    else
        optimalIndices(subj) = NaN; % Assign NaN if no peak is found in the interval
    end
    
    % Optional Plotting for selected subjects
    if ismember(subj, subjectsForPlot)
        % Create subplot for each selected subject
        subplot(ceil(length(subjectsForPlot) / 2), 2, find(subjectsForPlot == subj)); % Adjust for 2 columns layout
        plot(timePoints(:, 1), MEG_LI(subj, :), 'b'); % Plot MEG LI as a blue line
        hold on;
        
        % Highlight the optimal time point if it exists
        if ~isnan(optimalTimePointIdx)
            plot(timePoints(optimalTimePointIdx, 1), MEG_LI(subj, optimalTimePointIdx), 'ro'); % Mark optimal point as a red circle
        end
        
        % Set the y-axis limits
        ylim([-100 100]);
        
        hold off;
        if ~isnan(optimalTimePointIdx)
            title(sprintf('Subject %d | fMRI LI: %.2f | MEG LI: %.2f', subj, fMRI_LI(subj), MEG_LI(subj, optimalTimePointIdx)));
        else
            title(sprintf('Subject %d | fMRI LI: %.2f | No Optimal MEG LI Found', subj, fMRI_LI(subj)));
        end
        xlabel('Time Points');
        ylabel('MEG LI');
        set(gcf, 'Position', [200, 400, 800, 400]);
    end
end
end
