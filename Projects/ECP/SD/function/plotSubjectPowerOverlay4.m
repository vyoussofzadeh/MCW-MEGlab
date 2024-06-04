% function aucDiff = plotSubjectPowerOverlay4(powerLeft, powerRight, timePoints, lowerBound, upperBound)
% % Computes the area under the curve of the power difference between left and right hemispheres for all subjects within a specified interval.
%
% numSubjects = size(powerLeft, 1);
% aucDiff = zeros(numSubjects, 1); % Initialize array to store AUC for each subject
%
% % Create a figure to hold all subplots if plotting is enabled
% % if plot_flag == 1, figure; end
%
% % Identify the indices within the specified bounds
% timeIdx = find(timePoints >= lowerBound & timePoints <= upperBound);
%
% for subj = 1:numSubjects
%     % Calculate the difference in power between left and right hemispheres
%     powerDiff = powerLeft(subj, :) - powerRight(subj, :);
%
%     % Compute the area under the curve for the difference within bounds for each subject
%     aucDiff(subj) = trapz(timePoints(timeIdx), powerDiff(timeIdx));
%
% end
% end

function aucDiff = plotSubjectPowerOverlay4(powerLeft, powerRight, ID, timePoints, plot_flag, lowerBound, upperBound)
% Computes the area under the curve of the power difference between left and right hemispheres for all subjects within a specified interval.

numSubjects = size(powerLeft, 1);
aucDiff = zeros(numSubjects, 1); % Initialize array to store AUC for each subject

% Create a figure to hold all subplots if plotting is enabled
if plot_flag == 1
    figure; % Create a figure only if plotting is enabled
end

% Identify the indices within the specified bounds
timeIdx = find(timePoints >= lowerBound & timePoints <= upperBound);

for subj = 1:numSubjects
    % Calculate the difference in power between left and right hemispheres
    powerDiff = powerLeft(subj, :) - powerRight(subj, :);
    
    PL = powerLeft(subj, :);
    PR = powerRight(subj, :);
    
    % Compute the area under the curve for the difference within bounds for each subject
    aucDiff_L = trapz(timePoints(timeIdx), PL(timeIdx));
    aucDiff_R = trapz(timePoints(timeIdx), PR(timeIdx));
    aucDiff(subj) = aucDiff_L - aucDiff_R;
    
    if plot_flag == 1
        % Plotting each subject's power data and the difference
        subplot(ceil(numSubjects / 4), 4, subj); % Adjust layout based on the number of subjects
        hold on;
        plot(timePoints, powerLeft(subj, :), 'b', 'LineWidth', 1); % Plot left hemisphere power as a blue line
        plot(timePoints, powerRight(subj, :), 'r', 'LineWidth', 1); % Plot right hemisphere power as a red line
        
        % Create the fill vectors
        fillX = [timePoints(timeIdx), fliplr(timePoints(timeIdx))];
        fillY = [zeros(1, numel(timeIdx)), fliplr(powerDiff(timeIdx))];
        
        if length(fillX) == length(fillY)
            fill(fillX, fillY, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        else
            warning('fillX and fillY vectors are not of the same length');
        end
        
        hold off;
        axis tight;
        title(sprintf('Subject %s', ID{subj}));
        
        if subj == 1 % Add a legend only in the first subplot to avoid repetition
            legend({'Left Power', 'Right Power', 'Power Difference'}, 'Location', 'best');
        end
    end
end
xlabel('Time Points');
ylabel('Power (dB)');

if plot_flag == 1
    % Adjust figure properties
    set(gcf, 'Position', [100, 100, 1600, 900]); % Adjust figure size and position
end

end


% function aucDiff = plotSubjectPowerOverlay4(powerLeft, powerRight, ID, timePoints, plot_flag, lowerBound, upperBound)
% % Computes the area under the curve of the power difference between left and right hemispheres for all subjects within a specified interval.
%
% numSubjects = size(powerLeft, 1);
% aucDiff = zeros(numSubjects, 1); % Initialize array to store AUC for each subject
%
% % Create a figure to hold all subplots if plotting is enabled
% if plot_flag == 1
%     figure; % Create a figure only if plotting is enabled
% end
%
% % Identify the indices within the specified bounds
% timeIdx = find(timePoints >= lowerBound & timePoints <= upperBound);
%
% for subj = 1:numSubjects
%     % Calculate the difference in power between left and right hemispheres
%     powerDiff = powerLeft(subj, :) - powerRight(subj, :);
%
%     % Compute the area under the curve for the difference within bounds for each subject
%     aucDiff(subj) = trapz(timePoints(timeIdx), powerDiff(timeIdx));
%
%     if plot_flag == 1
%         % Plotting each subject's power data and the difference
%         subplot(ceil(numSubjects / 4), 4, subj); % Adjust layout based on the number of subjects
%         hold on;
%         plot(timePoints, powerLeft(subj, :), 'b', 'LineWidth', 1); % Plot left hemisphere power as a blue line
%         plot(timePoints, powerRight(subj, :), 'r', 'LineWidth', 1); % Plot right hemisphere power as a red line
%         fill([timePoints(timeIdx) fliplr(timePoints(timeIdx))], [zeros(size(powerDiff(timeIdx))) fliplr(powerDiff(timeIdx))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
%         hold off;
%         axis tight;
%         title(sprintf('Subject %s', ID{subj}));
%         xlabel('Time Points');
%         ylabel('Power (dB)');
%         if subj == 1 % Add a legend only in the first subplot to avoid repetition
%             legend({'Left Power', 'Right Power', 'Power Difference'}, 'Location', 'best');
%         end
%     end
% end
%
% if plot_flag == 1
%     % Adjust figure properties
%     set(gcf, 'Position', [100, 100, 1600, 900]); % Adjust figure size and position
% end
%
% end
