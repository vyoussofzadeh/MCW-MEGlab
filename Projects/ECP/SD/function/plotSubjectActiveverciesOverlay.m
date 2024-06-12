function [optimalTimePoints] = plotSubjectActiveverciesOverlay(roi_cnt_pt, ID, timePoints, plot_flag, lowerBound, upperBound)
% Plots the active vertices count for a region of interest (ROI) for all subjects with individual peak annotations.
% Searches for the optimal time point based on the maximum active vertices count strictly within the interval defined by lowerBound and upperBound.

numSubjects = size(roi_cnt_pt, 1);
optimalTimePoints = zeros(numSubjects, 1); % Initialize array to store optimal time points

% Create a figure to hold all subplots if plotting is enabled
if plot_flag == 1, figure; end

for subj = 1:numSubjects
    
    maxActiveVertices = 0; % Initialize maximum active vertices count
    optimalIdx = NaN; % Initialize index of optimal time point
    
    for i = 1:length(timePoints)
        % Check if the current time point is within the desired interval
        if timePoints(i) >= lowerBound && timePoints(i) <= upperBound
            % Find the time point with the maximum active vertices count within the interval
            if roi_cnt_pt(subj, i) > maxActiveVertices
                maxActiveVertices = roi_cnt_pt(subj, i);
                optimalIdx = i;
            end
        end
    end
    
    if ~isnan(optimalIdx)
        optimalTimePoints(subj) = timePoints(optimalIdx);
    else
        optimalTimePoints(subj) = NaN; % Assign NaN if no maximum is found in the interval
    end
    
    if plot_flag == 1
        subplot(ceil(numSubjects / 6), 6, subj); % Adjust layout to fit all subjects
        plot(timePoints, roi_cnt_pt(subj, :), 'LineWidth', 1.5); % Plot active vertices count
        
        % Highlight the optimal time point, if within bounds
        %         if optimalIdx > 0
        %             plot(optimalTimePoints(subj), roi_cnt_pt(subj, optimalIdx), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'y'); % Highlight the optimal time point
        %         end
        
        hold off;
        axis tight;
        title(sprintf('Subject %s', ID{subj}));
    end
end

if plot_flag == 1
    xlabel('Time Points');
    ylabel('Active Vertices Count');
    legend({'Active Vertices Count', 'Peak'}, 'Location', 'best', 'Orientation', 'horizontal');
    set(gcf, 'Position', [100   100   1200   1000]);
    
end

end


% function [optimalTimePoints] = plotSubjectActiveverciesOverlay(roi_cnt_pt, ID, timePoints, plot_flag, lowerBound, upperBound)
% % Plots overlaid power values for left and right hemispheres for all subjects with individual peak annotations.
% % Searches for the optimal time point based on the maximum absolute power difference strictly within the interval defined by lowerBound and upperBound.
%
% numSubjects = size(powerLeft, 1);
% optimalTimePoints = zeros(numSubjects, 1); % Initialize array to store optimal time points
%
% % Create a figure to hold all subplots if plotting is enabled
% if plot_flag == 1, figure; end
%
% for subj = 1:numSubjects
%
%     maxPowerDiff = 0; % Initialize maximum absolute MEG LI
%     optimalIdx = NaN; % Initialize index of optimal time point
%
%     %     % Calculate the difference in power between left and right hemispheres for plotting
%     powerDiff = powerLeft(subj, :) - powerRight(subj, :);
%
%     for i = 1:length(timePoints)
%         % Check if the current time point is within the desired interval
%         if timePoints(i) >= lowerBound && timePoints(i) <= upperBound
%             % Find the maximum absolute MEG LI within the interval
%             if (powerDiff(i)) > maxPowerDiff
%                 maxPowerDiff = (powerDiff(i));
%                 optimalIdx = i;
%             end
%         end
%     end
%
%
%     if ~isnan(optimalIdx)
%         optimalTimePoints(subj) = timePoints(optimalIdx);
%     else
%         optimalTimePoints(subj) = NaN; % Assign NaN if no maximum is found in the interval
%     end
%
%     if plot_flag == 1
%         subplot(ceil(numSubjects / 6), 6, subj); % Adjust layout to fit all subjects
%         plot(timePoints, powerLeft(subj, :), 'b'); % Plot left hemisphere power as a blue line
%         hold on;
%         plot(timePoints, powerRight(subj, :), 'r'); % Plot right hemisphere power as a red line
%         plot(timePoints, powerDiff, 'k--', 'LineWidth', 1.5); % Overlay the difference as a dashed black line
%
%         % Highlight the optimal time point, if within bounds
%         if optimalIdx > 0
%             plot(optimalTimePoints(subj), powerDiff(optimalIdx), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'y'); % Highlight the optimal time point
%         end
%
%         hold off;
%         axis tight;
%         title(sprintf('Subject %s', ID{subj}));
%     end
% end
%
% if plot_flag == 1
%     xlabel('Time Points');
%     ylabel('Power');
%     legend({'Left Power', 'Right Power', 'Left-Right Diff', 'Peak'}, 'Location', 'best', 'Orientation', 'horizontal');
% end
%
% end
%
%
