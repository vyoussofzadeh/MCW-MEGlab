function [optimalIndices, maxDiff_opt] = findIndividualOptimalTimePoints_interval_rSNR3(rSNR_left, rSNR_right, fMRI_LI, timePoints, subjectsForPlot, lowerBound, upperBound)

numSubjects = size(rSNR_left, 1);
optimalIndices = zeros(numSubjects, 1); % Initialize array for optimal indices
maxDiff_opt = zeros(numSubjects, 1); % Initialize array for maximum difference

% Create a figure outside the loop for all subplots
if ~isempty(subjectsForPlot) && mean(~isnan(subjectsForPlot))
    figure;
end

dff = ((rSNR_left(:)) - (rSNR_right(:)));
mx = max([max(dff(:)), max(rSNR_right(:)), max(rSNR_left(:))]);
mn = min([min(dff(:)), min(rSNR_right(:)), min(rSNR_left(:))]);

for subj = 1:numSubjects
    maxDiff = -inf; % Initialize maximum difference to negative infinity
    optimalTimePointIdx = NaN; % Initialize index of optimal time point
    timePointMid = (timePoints(:,1));
    
    for i = 1:length(timePointMid)
        % Calculate log squared ratio between rSNR_left and rSNR_right
        currentDiff = abs((rSNR_left(subj, i)) - (rSNR_right(subj, i)));
        
        
        % Check if the midpoint of the current time interval is within the desired interval
        if timePointMid(i) >= lowerBound && timePointMid(i) <= upperBound
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
            %             % Create subplot for each selected subject
            if length(subjectsForPlot) > 10
                subplot(ceil(length(subjectsForPlot) / 6), 6, find(subjectsForPlot == subj)); % Adjust for 6 columns layout
                set(gcf, 'Position', [100, 100, 1600, 1300]);
            else
                subplot(ceil(length(subjectsForPlot) / 2), 2, find(subjectsForPlot == subj)); % Adjust for 2 columns layout
                set(gcf, 'Position', [200, 400, 800, 400]);
            end
            p1 = plot(timePointMid, rSNR_left(subj, :), 'color', 'r', 'DisplayName', 'rSNR Left', 'LineWidth', 1.5); % Plot left SNR in Deep Teal
            hold on
            p2 = plot(timePointMid, rSNR_right(subj, :), 'color', 'b', 'DisplayName', 'rSNR Right', 'LineWidth', 1.5); % Plot right SNR in Warm Orange
            
%             yyaxis right; % Right Y-axis for log scale values
            set(gca, 'ycolor', 'k'); % Set y-axis color back to black for left axis
            val = ((rSNR_left(subj, :)) - (rSNR_right(subj, :)));
            
            p3 = plot(timePointMid, val, 'DisplayName', 'rSNR (L-R)', 'LineWidth', 1.5, 'color', [0.5,0.5,0.5]); % Plot difference as a blue line
            
            % Highlight the optimal time point
            plot(timePointMid(optimalTimePointIdx), val(optimalTimePointIdx), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimal Point'); % Mark optimal point as a black circle

            % Plot vertical lines for the lower and upper bounds
            xline(lowerBound, '--k', 'DisplayName', 'Lower Bound');
            xline(upperBound, '--k', 'DisplayName', 'Upper Bound');
            xlim([timePointMid(1) timePointMid(end)])
%             title(sprintf('S%d|Max Mean:%.1f', subj, maxDiff_opt(subj)));
            title(sprintf('S%d|fMRI:%.1f', subj, fMRI_LI(subj)));

            box off;
            set(gca, 'color', 'none');
%             ylim([mn, mx])

        end
    end
end

if ~isnan(subjectsForPlot)
    %     set(gcf, 'Position', [200, 400, 800, 400]); % Adjust size of figure window
    hold off;
    lgd = legend([p1, p2, p3], 'Location', 'best','Orientation', 'horizontal'); % Add a legend
    lgdPos = lgd.Position; % Get current position
    lgdPos(2) = lgdPos(2) - 0.08; % Move legend down
    lgdPos(1) = lgdPos(1) - 0.05; % Move legend down
    lgd.Position = lgdPos;
    ylabel('rSNR (L - R)');
    xlabel('Time Points');
    ylabel('rSNR (L & R)');

end

end



% function [optimalIndices, maxDiff_opt] = findIndividualOptimalTimePoints_interval_rSNR2(rSNR_left, rSNR_right, fMRI_LI, timePoints, subjectsForPlot, lowerBound, upperBound)
% 
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
% %     timePointMid = (timePoints(:,1) + timePoints(:,2)) / 2; % Midpoint of each time interval
%     timePointMid = (timePoints(:,1));
%     
%     for i = 1:length(timePointMid)
%         % Avoiding division by zero or extreme values with epsilon
%         %         epsilon = 1e-10;
%         % Calculate log squared ratio between rSNR_left and rSNR_right
%         %         currentDiff = log((rSNR_left(subj, i) + epsilon) / (rSNR_right(subj, i) + epsilon))^2;
%         %         currentDiff = ((rSNR_left(subj, i) + epsilon) / (rSNR_right(subj, i) + epsilon))^2;
%         currentDiff = abs((rSNR_left(subj, i)) - (rSNR_right(subj, i)));
%         
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
%             %             % Create subplot for each selected subject
%             if length(subjectsForPlot) > 10
%                 subplot(ceil(length(subjectsForPlot) / 6), 6, find(subjectsForPlot == subj)); % Adjust for 6 columns layout
%                 set(gcf, 'Position', [100, 100, 1600, 1300]);
%             else
%                 subplot(ceil(length(subjectsForPlot) / 2), 2, find(subjectsForPlot == subj)); % Adjust for 2 columns layout
%                 set(gcf, 'Position', [200, 400, 800, 400]);
%             end
%             yyaxis left; % Left Y-axis for original SNR values
%             p1 = plot(timePointMid, rSNR_left(subj, :), 'DisplayName', 'rSNR Left', 'LineWidth', 1.5); % Plot left SNR in Deep Teal
%             hold on
%             p2 = plot(timePointMid, rSNR_right(subj, :), 'color', 'r', 'DisplayName', 'rSNR Right', 'LineWidth', 1.5); % Plot right SNR in Warm Orange
%             
%             yyaxis right; % Right Y-axis for log scale values
%             set(gca, 'ycolor', 'k'); % Set y-axis color back to black for left axis
%             %             val = log((rSNR_left(subj, :) + epsilon) ./ (rSNR_right(subj, :) + epsilon)).^2;
%             %             val = ((rSNR_left(subj, :) + epsilon) ./ (rSNR_right(subj, :) + epsilon)).^2;
%             val = ((rSNR_left(subj, :)) - (rSNR_right(subj, :)));
%             
%             p3 = plot(timePointMid, val, 'DisplayName', 'rSNR (L-R)', 'LineWidth', 1.5, 'color', 'k'); % Plot difference as a blue line
%             
%             % Highlight the optimal time point
%             scatter(timePointMid(optimalTimePointIdx), val(optimalTimePointIdx), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimal Point'); % Mark optimal point as a black circle
%             
%             % Plot vertical lines for the lower and upper bounds
%             xline(lowerBound, '--k', 'DisplayName', 'Lower Bound');
%             xline(upperBound, '--k', 'DisplayName', 'Upper Bound');
%             xlim([timePointMid(1) timePointMid(end)])
%             yyaxis left; % Switch back to left axis for common settings
% %             title(sprintf('S%d|Max Mean:%.1f', subj, maxDiff_opt(subj)));
%             title(sprintf('S%d|fMRI:%.1f', subj, fMRI_LI(subj)));
% 
%             box off;
%             set(gca, 'color', 'none');
%         end
%     end
% end
% 
% if ~isnan(subjectsForPlot)
%     %     set(gcf, 'Position', [200, 400, 800, 400]); % Adjust size of figure window
%     hold off;
%     lgd = legend([p1, p2, p3], 'Location', 'best','Orientation', 'horizontal'); % Add a legend
%     lgdPos = lgd.Position; % Get current position
%     lgdPos(2) = lgdPos(2) - 0.08; % Move legend down
%     lgdPos(1) = lgdPos(1) - 0.05; % Move legend down
%     %     lgdPos(1) = lgdPos(1) + 0.05; % Move legend down
%     lgd.Position = lgdPos;
%     yyaxis right; % Left Y-axis for original SNR values
%     ylabel('rSNR (L - R)');
%     xlabel('Time Points');
%     yyaxis left; % Left Y-axis for original SNR values
%     ylabel('rSNR (L & R)');
% 
% end
% 
% end
% 
