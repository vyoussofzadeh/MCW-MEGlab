function [optimalIndices, maxDiff_opt, rsnr_max] = findIndividualOptimalTimePoints_interval_rSNR4(rSNR_left, rSNR_right, fMRI_LI, timePoints, subjectsForPlot, bounds)

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
    timePointMid = (timePoints(:,1));
    
    rsnr_max_within_bounds = -inf; % Initialize the max rSNR within bounds to negative infinity

    for i = 1:length(timePointMid)
        % Calculate log squared ratio between rSNR_left and rSNR_right
        currentDiff = abs((rSNR_left(subj, i)) - (rSNR_right(subj, i)));
        
        % Check if the midpoint of the current time interval is within the desired interval
        if timePointMid(i) >= timePointMid(bounds(subj,1)) && timePointMid(i) <= timePointMid(bounds(subj,2))
            % Find the maximum difference within the interval
            if currentDiff > maxDiff
                maxDiff = currentDiff;
                optimalTimePointIdx = i;
                maxDiff_opt(subj) = currentDiff;
            end
            % Update the maximum rSNR within bounds
            rsnr_max_within_bounds = max(rsnr_max_within_bounds, max(rSNR_left(subj, i), rSNR_right(subj, i)));
        end
    end
    
    if ~isnan(optimalTimePointIdx)
        optimalIndices(subj) = optimalTimePointIdx; % Store the index of the optimal time point
    else
        optimalIndices(subj) = NaN; % Assign NaN if no maximum is found in the interval
    end
    
    rsnr_max(subj) = rsnr_max_within_bounds; % Assign the maximum rSNR found within bounds

    if ~isnan(subjectsForPlot)
        % Optional Plotting for selected subjects
        if ismember(subj, subjectsForPlot)
            % Create subplot for each selected subject
            if length(subjectsForPlot) > 10
                subplot(ceil(length(subjectsForPlot) / 6), 6, find(subjectsForPlot == subj)); % Adjust for 6 columns layout
                set(gcf, 'Position', [100, 100, 1600, 1300]);
            else
                subplot(ceil(length(subjectsForPlot) / 2), 2, find(subjectsForPlot == subj)); % Adjust for 2 columns layout
                set(gcf, 'Position', [200, 400, 800, 400]);
            end
            plot(timePointMid, rSNR_left(subj, :), 'color', 'r', 'DisplayName', 'rSNR Left', 'LineWidth', 1.5);
            hold on
            plot(timePointMid, rSNR_right(subj, :), 'color', 'b', 'DisplayName', 'rSNR Right', 'LineWidth', 1.5);
            plot(timePointMid, rSNR_left(subj, :) - rSNR_right(subj, :), 'color', [0.5,0.5,0.5], 'DisplayName', 'rSNR (L-R)', 'LineWidth', 1.5);
            
            % Highlight the optimal time point
            plot(timePointMid(optimalTimePointIdx), rSNR_left(subj, optimalTimePointIdx) - rSNR_right(subj, optimalTimePointIdx), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimal Point'); % Mark optimal point as a black circle
            
            % Plot vertical lines for the lower and upper bounds
            xline(timePointMid(bounds(subj,1)), '--k', 'DisplayName', 'Lower Bound');
            xline(timePointMid(bounds(subj,2)), '--k', 'DisplayName', 'Upper Bound');
            xlim([timePointMid(1) timePointMid(end)])
            title(sprintf('S%d|fMRI:%.1f', subj, fMRI_LI(subj)));
            
            box off;
            set(gca, 'color', 'none');
        end
    end
end

if ~isnan(subjectsForPlot)
    hold off;
    legend('Location', 'best', 'Orientation', 'horizontal');
    xlabel('Time Points');
    ylabel('rSNR (L & R)');
end
end

