function [optimalIndices, maxDiff_opt, bounds] = findIndividualOptimalTimePoints_interval_rSNR_optbound(rSNR_left, rSNR_right, timePoints, subjectsForPlot, minlowerband, maxUpperband)

% Identify non-negative timePoints indices
positiveIndices = find(timePoints(:,1) >= minlowerband & timePoints(:,1) <= maxUpperband);

filteredTimePoints = timePoints(positiveIndices,:);
rSNR_left_pos = rSNR_left(:, positiveIndices);
rSNR_right_pos = rSNR_right(:, positiveIndices);

% Initialize variables
numSubjects = size(rSNR_left_pos, 1);
optimalIndices = zeros(numSubjects, 1);
maxDiff_opt = zeros(numSubjects, 1);
bounds = zeros(numSubjects, 2);

% Calculate absolute differences for both original and filtered data
abs_diff = abs(rSNR_left - rSNR_right);
abs_diff_pos = abs(rSNR_left_pos - rSNR_right_pos);

% Calculate percentiles and thresholds for filtered data
percentile95 = prctile(abs_diff_pos, 95, 2);
thresholds = 0.05 * max(percentile95, [], 2);

% Create a figure outside the loop for all subplots
if ~isempty(subjectsForPlot) && mean(~isnan(subjectsForPlot))
    figure;
end

for subj = 1:numSubjects
    maxDiff = -inf;
    optimalTimePointIdx = NaN;
    firstTimePointIdx = NaN;
    lastTimePointIdx = NaN;

    for i = 1:length(filteredTimePoints)
        currentDiff = abs_diff(subj, positiveIndices(i));

        % Check if the current difference exceeds the dynamic threshold
        if currentDiff > thresholds(subj)
            if isnan(firstTimePointIdx)
                firstTimePointIdx = i;
            end
            lastTimePointIdx = i;

            if currentDiff > maxDiff
                maxDiff = currentDiff;
                optimalTimePointIdx = i;
                maxDiff_opt(subj) = currentDiff;
            end
        end
    end

    if ~isnan(optimalTimePointIdx)
        optimalIndices(subj) = positiveIndices(optimalTimePointIdx);
        bounds(subj, :) = positiveIndices([firstTimePointIdx, lastTimePointIdx]);
    else
        optimalIndices(subj) = NaN;
        bounds(subj, :) = [NaN, NaN];
    end

    % Optional plotting for selected subjects
    if ~isempty(subjectsForPlot) && ismember(subj, subjectsForPlot)
        plotSubjectData(subj, subjectsForPlot, timePoints(:,1), rSNR_left(subj, :), rSNR_right(subj, :), abs_diff(subj, :), optimalIndices(subj), timePoints(bounds(subj, 1), 1), timePoints(bounds(subj, 2), 1), thresholds(subj));
    end
end

end

function plotSubjectData(subj, subjectsForPlot, timePointMid, rSNR_left, rSNR_right, abs_diff, optimalIdx, lowerBound, upperBound, thresholds)
subplot(ceil(length(subjectsForPlot) / 2), 2, find(subjectsForPlot == subj));
set(gcf, 'Position', [200, 400, 800, 400]);
%     yyaxis left;
plot(timePointMid, rSNR_left, 'b-', 'DisplayName', 'rSNR Left');
hold on;
plot(timePointMid, rSNR_right, 'r-', 'DisplayName', 'rSNR Right');
xline(lowerBound, '--k', 'DisplayName', 'Lower Bound');
xline(upperBound, '--k', 'DisplayName', 'Upper Bound');

yline(thresholds)

%     yyaxis right;
plot(timePointMid, abs_diff, 'k-', 'DisplayName', 'Absolute Difference', 'LineWidth', 2);
scatter(timePointMid(optimalIdx), abs_diff(optimalIdx), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimal Point');

title(sprintf('Subject %d | Max Difference: %.1f', subj, abs_diff(optimalIdx)));
%     legend('show');
box off;
set(gca, 'color', 'none');
end
