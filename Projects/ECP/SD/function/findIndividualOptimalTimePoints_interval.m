function [optimalIndices, maxLI_opt] = findIndividualOptimalTimePoints_interval(MEG_LI, fMRI_LI, timePoints, subjectsForPlot, lowerBound, upperBound)
% This function finds the optimal time points based on the maximum absolute MEG laterality index within specified bounds.

% Inputs:
%   MEG_LI - Matrix of MEG laterality indices for each subject over time
%   fMRI_LI - Laterality indices from fMRI for comparison or additional context
%   timePoints - Matrix where each row contains the start and end of each time interval
%   subjectsForPlot - Indices of subjects for whom plots should be generated
%   lowerBound, upperBound - Time bounds within which to find optimal time points

% Outputs:
%   optimalIndices - Indices of the optimal time points for each subject
%   maxLI_opt - Maximum absolute LI values at these optimal time points

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
        % Compute the absolute MEG LI for the current time point
        current_LI = abs(MEG_LI(subj, i));
        
        % Check if the midpoint of the current time interval is within the desired interval
        if timePointMid(i) >= lowerBound && timePointMid(i) <= upperBound
            % Find the maximum absolute LI within the interval
            if current_LI > maxMEG_LI
                maxMEG_LI = current_LI;
                optimalTimePointIdx = i;
                maxLI_opt(subj) = current_LI;
            end
        end
    end
    
    % Store the index of the optimal time point
    if ~isnan(optimalTimePointIdx)
        optimalIndices(subj) = optimalTimePointIdx;
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
            plot(timePointMid, (MEG_LI(subj, :)), 'b'); % Plot absolute MEG LI as a blue line
            hold on;
            
            % Highlight the optimal time point
            plot(timePointMid(optimalTimePointIdx), maxLI_opt(subj), 'ro','MarkerSize', 6, 'MarkerFaceColor', 'r'); % Mark optimal point as a red circle
            
            % Set the y-axis limits
            ylim([-100 100]); % Adjusted to only positive values as we are dealing with absolute values
            
            % Plot vertical lines for the lower and upper bounds
            xline(lowerBound, '--k');
            xline(upperBound, '--k');
            
            hold off;
            set(gca, 'color', 'none');
            title(sprintf('S%d|fMRI:%.1f|MEG:%.1f', subj, fMRI_LI(subj), maxLI_opt(subj)));
            
            set(gcf, 'Position', [200, 400, 800, 400]);
            box off
        end
    end
end

if ~isnan(subjectsForPlot)
    xlabel('Time Points');
    ylabel('Abs(MEG LI)');
end

end


% function [optimalIndices, maxLI_opt] = findIndividualOptimalTimePoints_interval(MEG_LI, fMRI_LI, timePoints, subjectsForPlot, lowerBound, upperBound)
% 
% numSubjects = size(MEG_LI, 1);
% optimalIndices = zeros(numSubjects, 1); % Initialize array for optimal indices
% maxLI_opt = zeros(numSubjects, 1); % Initialize array for maximum LI
% 
% % Create a figure outside the loop for all subplots
% if ~isempty(subjectsForPlot) && mean(~isnan(subjectsForPlot))
%     figure;
% end
% 
% for subj = 1:numSubjects
%     maxMEG_LI = -inf; % Initialize maximum MEG LI to negative infinity
%     optimalTimePointIdx = NaN; % Initialize index of optimal time point
%     timePointMid = (timePoints(:,1) + timePoints(:,2)) / 2; % Midpoint of each time interval
%     
%     for i = 1:length(timePointMid)
%         % Check if the midpoint of the current time interval is within the desired interval
%         if timePointMid(i) >= lowerBound && timePointMid(i) <= upperBound
%             % Find the maximum MEG LI within the interval
%             if MEG_LI(subj, i) > maxMEG_LI
%                 maxMEG_LI = MEG_LI(subj, i);
%                 optimalTimePointIdx = i;
%                 maxLI_opt(subj) = MEG_LI(subj, i);
%             end
%         end
%     end
%     
%     if ~isnan(optimalTimePointIdx)
%         optimalIndices(subj) = optimalTimePointIdx; % Store the index of the optimal time point
%     else
%         optimalIndices(subj) = NaN; % Assign NaN if no maximum is found in the interval
%     end
%     
%     if ~isnan(subjectsForPlot)
%         % Optional Plotting for selected subjects
%         if ismember(subj, subjectsForPlot)
%             % Create subplot for each selected subject
%             if length(subjectsForPlot) > 20
%                 subplot(ceil(length(subjectsForPlot) / 6), 6, find(subjectsForPlot == subj)); % Adjust for 6 columns layout
%             else
%                 subplot(ceil(length(subjectsForPlot) / 2), 2, find(subjectsForPlot == subj)); % Adjust for 2 columns layout
%             end
%             plot(timePointMid, MEG_LI(subj, :), 'b'); % Plot MEG LI as a blue line
%             hold on;
%             
%             % Highlight the optimal time point
%             plot(timePointMid(optimalTimePointIdx), MEG_LI(subj, optimalTimePointIdx), 'ro','MarkerSize', 6, 'MarkerFaceColor', 'r'); % Mark optimal point as a red circle
%             
%             
%             % Set the y-axis limits
%             ylim([-100 100]);
%             
%             % Plot vertical lines for the lower and upper bounds
%             xline(lowerBound, '--k');
%             xline(upperBound, '--k');
%             
%             hold off;
%             set(gca, 'color', 'none');
%             %         title(sprintf('S%d|fMRI LI:%.2f|MEG LI:%.2f', subj, fMRI_LI(subj), MEG_LI(subj, optimalTimePointIdx)));
%             title(sprintf('S%d|fMRI:%.1f|MEG:%.1f', subj, fMRI_LI(subj), MEG_LI(subj, optimalTimePointIdx)));
%             
%             set(gcf, 'Position', [200, 400, 800, 400]);
%             box off
%         end
%     end
% end
% 
% if ~isnan(subjectsForPlot)
%     xlabel('Time Points');
%     ylabel('MEG LI');
% end
% % Set common xlabel and ylabel for the entire figure
% % han = gcf;
% % if ~isempty(han.Children)
% %     han.Children(end).XLabel.String = 'Time Points';
% %     han.Children(end).YLabel.String = 'MEG LI';
% % end
% 
% end


% function [optimalIndices, maxLI_opt] = findIndividualOptimalTimePoints_interval(MEG_LI, fMRI_LI, timePoints, subjectsForPlot, lowerBound, upperBound)
% 
% numSubjects = size(MEG_LI, 1);
% optimalIndices = NaN(numSubjects, 1); % Initialize array for optimal indices with NaN
% maxLI_opt = NaN(numSubjects, 1); % Initialize array for maximum LI with NaN
% 
% % Create a figure outside the loop for all subplots
% if ~isempty(subjectsForPlot) && any(~isnan(subjectsForPlot))
%     figure;
%     % Set the figure position and size more universally
%     set(gcf, 'Position', [100, 100, 1200, 800]);
% end
% 
% for subj = 1:numSubjects
%     maxMEG_LI = -inf; % Initialize maximum MEG LI to negative infinity
%     optimalTimePointIdx = NaN; % Initialize index of optimal time point
%     timePointMid = (timePoints(:,1) + timePoints(:,2)) / 2; % Midpoint of each time interval
% 
%     for i = 1:length(timePointMid)
%         % Check if the midpoint of the current time interval is within the desired interval
%         if timePointMid(i) >= lowerBound && timePointMid(i) <= upperBound
%             % Find the maximum MEG LI within the interval
%             if MEG_LI(subj, i) > maxMEG_LI
%                 maxMEG_LI = MEG_LI(subj, i);
%                 optimalTimePointIdx = i;
%                 maxLI_opt(subj) = MEG_LI(subj, i);
%             end
%         end
%     end
% 
%     if ~isnan(optimalTimePointIdx)
%         optimalIndices(subj) = optimalTimePointIdx; % Store the index of the optimal time point
%     end
% 
%     if ~isnan(subjectsForPlot) && ismember(subj, subjectsForPlot)
%         % Optional Plotting for selected subjects
%         subplotIdx = find(subjectsForPlot == subj);
%         if length(subjectsForPlot) > 20
%             subplot(ceil(length(subjectsForPlot) / 6), 6, subplotIdx);
%         else
%             subplot(ceil(length(subjectsForPlot) / 2), 2, subplotIdx);
%         end
%         plot(timePointMid, MEG_LI(subj, :), 'b'); % Plot MEG LI as a blue line
%         hold on;
%         
%         % Highlight the optimal time point
%         scatter(timePointMid(optimalTimePointIdx), MEG_LI(subj, optimalTimePointIdx), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); % Mark optimal point as a red circle
%         
%         % Plot vertical lines for the lower and upper bounds
%         xline(lowerBound, '--k', 'Label', 'Lower Bound');
%         xline(upperBound, '--k', 'Label', 'Upper Bound');
% 
%         hold off;
%         set(gca, 'color', 'none');
%         ylim([-100 100]); % Set the y-axis limits
%         title(sprintf('S%d | fMRI: %.1f | MEG: %.1f', subj, fMRI_LI(subj), maxLI_opt(subj)));
%         box off;
%     end
% end
% 
% if ~isempty(subjectsForPlot)
%     % Setting common xlabel and ylabel for the entire figure after plotting is done
%     supxlabel('Time Points');
%     supylabel('MEG LI');
% end
% 
% end


