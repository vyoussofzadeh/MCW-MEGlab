% function [optimalIndices, maxDiff_opt] = findIndividualOptimalTimePoints_interval_rSNR_MaxInd(rSNR_left, rSNR_right, timePoints, subjectsForPlot, lowerBound, upperBound)


function [optimalIndices, maxDiff_opt] = findIndividualOptimalTimePoints_interval_rSNR_MaxInd(rSNR_left, rSNR_right, timePoints, subjectsForPlot, lowerBound, upperBound)

numSubjects = size(rSNR_left, 1);
optimalIndices = zeros(numSubjects, 1); % Initialize array for optimal indices
maxDiff_opt = zeros(numSubjects, 1); % Initialize array for maximum difference

% Create a figure outside the loop for all subplots
if ~isempty(subjectsForPlot) && mean(~isnan(subjectsForPlot))
    figure;
end

for subj = 1:numSubjects
    maxDiff = -inf; % Initialize maximum difference to negative infinity
    optimalTimePointIdx = NaN; % Initialize index of optimal time point
    timePointMid = (timePoints(:,1) + timePoints(:,2)) / 2; % Midpoint of each time interval

    % Normalize data
    normalizedLeft = (rSNR_left(subj, :) - mean(rSNR_left(subj, :))) / std(rSNR_left(subj, :));
    normalizedRight = (rSNR_right(subj, :) - mean(rSNR_right(subj, :))) / std(rSNR_right(subj, :));

    for i = 1:length(timePointMid)
        % Calculate absolute difference of normalized data
        if timePointMid(i) >= lowerBound && timePointMid(i) <= upperBound
            currentDiff = abs(normalizedLeft(i) - normalizedRight(i));
            
            % Find the maximum difference within the interval
            if currentDiff > maxDiff
                maxDiff = currentDiff;
                optimalTimePointIdx = i;
                maxDiff_opt(subj) = currentDiff;
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
            subplot(ceil(length(subjectsForPlot) / 6), 6, find(subjectsForPlot == subj)); % Adjust for 6 columns layout
%             subplot(subplotLayout(1), subplotLayout(2), find(subjectsForPlot == subj)); % Adjust for the layout
            
            hold on;
            % Plot normalized SNR lines
            plot(timePointMid, normalizedLeft, 'g-', 'DisplayName', 'Normalized rSNR Left'); % Plot normalized left SNR in green
            plot(timePointMid, normalizedRight, 'r-', 'DisplayName', 'Normalized rSNR Right'); % Plot normalized right SNR in red
            
            % Highlight the optimal time point
            scatter(timePointMid(optimalTimePointIdx), maxDiff_opt(subj), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimal Point'); % Mark optimal point as a black circle
            
            % Plot vertical lines for the lower and upper bounds
            xline(lowerBound, '--k', 'DisplayName', 'Lower Bound');
            xline(upperBound, '--k', 'DisplayName', 'Upper Bound');
            
            hold off;
            title(sprintf('S%d|Max Diff:%.2f', subj, maxDiff_opt(subj)));
            set(gca, 'color', 'none');
            legend('show');
            xlabel('Time Points');
            ylabel('Normalized rSNR');
            box off;
            set(gcf, 'Position', [200, 400, 800, 400]); % Adjust size of figure window
        end
    end
end

if ~isnan(subjectsForPlot)
    xlabel('Time Points');
    ylabel('Normalized rSNR Difference');
end

end

% numSubjects = size(rSNR_left, 1);
% optimalIndices = zeros(numSubjects, 1); % Initialize array for optimal indices
% maxDiff_opt = zeros(numSubjects, 1); % Initialize array for maximum difference
% 
% % Create a figure outside the loop for all subplots
% if ~isempty(subjectsForPlot) && mean(~isnan(subjectsForPlot))
%     figure;
% end
% 
% for subj = 1:numSubjects
%     maxDiff = -inf; % Initialize maximum difference to negative infinity
%     optimalTimePointIdx = NaN; % Initialize index of optimal time point
%     timePointMid = (timePoints(:,1) + timePoints(:,2)) / 2; % Midpoint of each time interval
%     
%     for i = 1:length(timePointMid)
%         % Calculate log squared ratio between rSNR_left and rSNR_right
%         % Avoid division by zero or negative values by adding a small constant if necessary
%         if rSNR_right(subj, i) ~= 0
%             currentDiff = log(rSNR_left(subj, i) / rSNR_right(subj, i))^2;
%         else
%             currentDiff = log((rSNR_left(subj, i) + 1e-10) / (rSNR_right(subj, i) + 1e-10))^2;
%         end
%         
%         % Check if the midpoint of the current time interval is within the desired interval
%         if timePointMid(i) >= lowerBound && timePointMid(i) <= upperBound
%             % Find the maximum difference within the interval
%             if currentDiff > maxDiff
%                 maxDiff = currentDiff;
%                 optimalTimePointIdx = i;
%                 maxDiff_opt(subj) = currentDiff;
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
%                 set(gcf, 'Position', [100, 100, 1200, 1000]);
%             else
%                 subplot(ceil(length(subjectsForPlot) / 2), 2, find(subjectsForPlot == subj)); % Adjust for 2 columns layout
%                 set(gcf, 'Position', [200, 400, 800, 400]);
%             end
%             
%             hold on;
%             % Plot individual SNR lines
%             plot(timePointMid, rSNR_left(subj, :), 'g-', 'DisplayName', 'rSNR Left'); % Plot left SNR in green
%             plot(timePointMid, rSNR_right(subj, :), 'r-', 'DisplayName', 'rSNR Right'); % Plot right SNR in red
%             plot(timePointMid, log(rSNR_left(subj, :) ./ rSNR_right(subj, :) + 1e-10).^2, 'b', 'DisplayName', 'Log square ratio (Left / Right)'); % Plot log squared ratio as a blue line
%             
%             % Highlight the optimal time point
%             scatter(timePointMid(optimalTimePointIdx), maxDiff_opt(subj), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimal Point'); % Mark optimal point as a black circle
%             
%             % Plot vertical lines for the lower and upper bounds
%             xline(lowerBound, '--k', 'DisplayName', 'Lower Bound');
%             xline(upperBound, '--k', 'DisplayName', 'Upper Bound');
%             
%             hold off;
%             title(sprintf('S%d|Max Mean:%.1f', subj, maxDiff_opt(subj)));
%             box off
%         end
%     end
% end
% 
% if ~isnan(subjectsForPlot)
%     legend('show', 'Location', 'best', 'Orientation', 'horizontal'); % Add a legend
%     xlabel('Time Points');
%     ylabel('rSNR');
% end
% 
% end
