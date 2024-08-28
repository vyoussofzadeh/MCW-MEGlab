function [optimalIndices, maxLI_opt] = findIndividualOptimalTimePoints_interval_edit(MEG_LI, fMRI_LI, timePoints, subjectsForPlot, lowerBound, upperBound)

numSubjects = size(MEG_LI, 1);
optimalIndices = zeros(numSubjects, 1); % Initialize array for optimal indices
maxLI_opt = zeros(numSubjects, 1); % Initialize array for maximum LI

% Create a figure outside the loop for all subplots
if ~isempty(subjectsForPlot) && mean(~isnan(subjectsForPlot))
    figure;
end

for subj = 1:numSubjects
    maxMEG_LI = -inf; % Initialize maximum MEG LI to negative infinity
    optimalTimePointIdx = NaN; % Initialize index of optimal time point
    timePointMid = (timePoints(:,1) + timePoints(:,2)) / 2; % Midpoint of each time interval
    
    for i = 1:length(timePointMid)
        % Check if the midpoint of the current time interval is within the desired interval
        if timePointMid(i) >= lowerBound && timePointMid(i) <= upperBound
            % Find the maximum MEG LI within the interval
            if MEG_LI(subj, i) > maxMEG_LI
                maxMEG_LI = MEG_LI(subj, i);
                optimalTimePointIdx = i;
                maxLI_opt(subj) = MEG_LI(subj, i);
            end
        end
    end
    
    if ~isnan(optimalTimePointIdx)
        optimalIndices(subj) = optimalTimePointIdx; % Store the index of the optimal time point
    else
        optimalIndices(subj) = NaN; % Assign NaN if no maximum is found in the interval
    end
    
    if ~isnan(subjectsForPlot)
        % Optional Plotting for selected subjects
        if ismember(subj, subjectsForPlot)
            % Create subplot for each selected subject
            if length(subjectsForPlot) > 20
                subplot(ceil(length(subjectsForPlot) / 6), 6, find(subjectsForPlot == subj)); % Adjust for 6 columns layout
            else
                subplot(ceil(length(subjectsForPlot) / 2), 2, find(subjectsForPlot == subj)); % Adjust for 2 columns layout
            end
            plot(timePointMid, MEG_LI(subj, :), 'b'); % Plot MEG LI as a blue line
            hold on;
            
            % Highlight the optimal time point
            plot(timePointMid(optimalTimePointIdx), MEG_LI(subj, optimalTimePointIdx), 'ro'); % Mark optimal point as a red circle
            
            
            % Set the y-axis limits
            ylim([-100 100]);
            
            % Plot vertical lines for the lower and upper bounds
            xline(lowerBound, '--k');
            xline(upperBound, '--k');
            
            hold off;
            set(gca, 'color', 'none');
            %         title(sprintf('S%d|fMRI LI:%.2f|MEG LI:%.2f', subj, fMRI_LI(subj), MEG_LI(subj, optimalTimePointIdx)));
            title(sprintf('S%d|fMRI:%.1f|MEG:%.1f', subj, fMRI_LI(subj), MEG_LI(subj, optimalTimePointIdx)));
            
            set(gcf, 'Position', [200, 400, 800, 400]);
        end
    end
end

if ~isnan(subjectsForPlot)
    xlabel('Time Points');
    ylabel('MEG LI');
end
% Set common xlabel and ylabel for the entire figure
% han = gcf;
% if ~isempty(han.Children)
%     han.Children(end).XLabel.String = 'Time Points';
%     han.Children(end).YLabel.String = 'MEG LI';
% end

end
