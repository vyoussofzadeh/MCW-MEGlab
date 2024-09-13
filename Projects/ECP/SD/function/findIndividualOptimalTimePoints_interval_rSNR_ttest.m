function [optimalIndices, minPValue] = findIndividualOptimalTimePoints_interval_rSNR_ttest(rSNR_left, rSNR_right, timePoints, subjectsForPlot, lowerBound, upperBound)

numSubjects = size(rSNR_left, 1);
optimalIndices = NaN(numSubjects, 1); % Initialize array for optimal indices with NaN
minPValue = ones(numSubjects, 1); % Initialize array for minimum p-values, starting at 1 (no significance)

% Create a figure outside the loop for all subplots
if ~isempty(subjectsForPlot) && mean(~isnan(subjectsForPlot))
    figure;
end

for subj = 1:numSubjects
    % Find indices within the time bounds that also do not exceed the length of the SNR arrays
    idxWithinBounds = find(timePoints >= lowerBound & timePoints <= upperBound & timePoints <= length(rSNR_left(subj, :)));

    if ~isempty(idxWithinBounds)
        % Ensure rSNR arrays for the subject are not shorter than the max index in idxWithinBounds
        if max(idxWithinBounds) <= length(rSNR_left(subj, :))
            % Conduct t-test within the bounds
            [h, pValues] = ttest(rSNR_left(subj, idxWithinBounds), rSNR_right(subj, idxWithinBounds));

            % Check each hypothesis test result and update the optimal index where the smallest significant p-value occurs
            significantIndices = find(h);
            if ~isempty(significantIndices)
                [minP, idx] = min(pValues(significantIndices));
                if minP < minPValue(subj)
                    minPValue(subj) = minP;
                    optimalIndices(subj) = idxWithinBounds(significantIndices(idx)); % Map back to the original array's index
                end
            end
        end
    end

    if ~isempty(subjectsForPlot) && ismember(subj, subjectsForPlot)
        % Optional Plotting for selected subjects
        subplotLayout = ceil(length(subjectsForPlot) / 6);
        subplot(subplotLayout, 6, find(subjectsForPlot == subj)); % Adjust for 6 columns layout
        
        hold on;
        % Plot rSNR lines
        plot(timePoints, rSNR_left(subj, :), 'g-', 'DisplayName', 'rSNR Left'); % Plot left SNR in green
        plot(timePoints, rSNR_right(subj, :), 'r-', 'DisplayName', 'rSNR Right'); % Plot right SNR in red
        
        if ~isnan(optimalIndices(subj))
            % Highlight the optimal time point window
            scatter(timePoints(optimalIndices(subj)), 0, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimal Point'); % Mark optimal point as a black circle
            title(sprintf('S%d|Min P-Value:%.4f', subj, minPValue(subj)));
        end
        
        hold off;
        legend('show');
        xlabel('Time Points');
        ylabel('rSNR');
        box off;
        set(gcf, 'Position', [200, 400, 800, 400]); % Adjust size of figure window
    end
end

if ~isnan(subjectsForPlot)
    xlabel('Time Points');
    ylabel('rSNR');
end

end
