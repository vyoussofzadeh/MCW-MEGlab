function [optimalIndices, maxAUC_opt] = findIndividualOptimalTimePoints_maxAUC(rSNR_left, rSNR_right, timePoints, subjectsForPlot, lowerBound, upperBound)

numSubjects = size(rSNR_left, 1);
optimalIndices = zeros(numSubjects, 1); % Initialize array for optimal indices
maxAUC_opt = zeros(numSubjects, 1); % Initialize array for maximum AUC

% Create a figure outside the loop for all subplots
if ~isempty(subjectsForPlot) && mean(~isnan(subjectsForPlot))
    figure('Position', [100, 100, 1200, 800]); % Adjust the figure size for better visibility
end

for subj = 1:numSubjects
    maxAUC = -inf; % Initialize maximum AUC to negative infinity
%     optimalTimePointIdx = NaN; % Initialize index of optimal time point
    timePointMid = (timePoints(:,1) + timePoints(:,2)) / 2; % Midpoint of each time interval
    
    % Determine the indices within the desired time bounds
    idxWithinBounds = find(timePointMid >= lowerBound & timePointMid <= upperBound);
    
    if ~isempty(idxWithinBounds)
        % Calculate the AUC for the maximum SNR between left and right within the bounds
        maxSNR = max(rSNR_left(subj, idxWithinBounds), rSNR_right(subj, idxWithinBounds));
        currentAUC = trapz(timePointMid(idxWithinBounds), maxSNR);
        
        % Compare with the maximum found so far
        if currentAUC > maxAUC
            maxAUC = currentAUC;
            optimalIndices(subj) = idxWithinBounds(1); % Store the index of the start of the optimal time point range
            maxAUC_opt(subj) = maxAUC;
        end
    end
    
    if ~isempty(subjectsForPlot) && ismember(subj, subjectsForPlot)
        % Optional Plotting for selected subjects
        subplotSetup = ceil(length(subjectsForPlot) / 6);
        subplot(subplotSetup, 6, find(subjectsForPlot == subj)); % Adjust for 6 columns layout

        hold on;
        % Plot individual SNR lines
        plot(timePointMid, rSNR_left(subj, :), 'g-', 'DisplayName', 'rSNR Left'); % Plot left SNR in green
        plot(timePointMid, rSNR_right(subj, :), 'r-', 'DisplayName', 'rSNR Right'); % Plot right SNR in red

        if ~isnan(optimalIndices(subj))
            % Highlight the optimal time point range
            scatter(timePointMid(optimalIndices(subj)), maxAUC_opt(subj), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimal AUC Start'); % Mark the start of optimal AUC as a black circle
        end

        % Plot vertical lines for the lower and upper bounds
        xline(lowerBound, '--k', 'DisplayName', 'Lower Bound');
        xline(upperBound, '--k', 'DisplayName', 'Upper Bound');

        hold off;
        title(sprintf('S%d|Max AUC:%.1f', subj, maxAUC_opt(subj)));
        legend('show', 'Location', 'best');
        xlabel('Time Points');
        ylabel('Max SNR');
    end
end
end
