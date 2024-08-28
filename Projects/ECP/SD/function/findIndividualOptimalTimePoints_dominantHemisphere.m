function [optimalIndices, maxSNR_opt] = findIndividualOptimalTimePoints_dominantHemisphere(rSNR_left, rSNR_right, timePoints, subjectsForPlot, lowerBound, upperBound)

numSubjects = size(rSNR_left, 1);
optimalIndices = zeros(numSubjects, 1); % Initialize array for optimal indices
maxSNR_opt = zeros(numSubjects, 1); % Initialize array for maximum SNR

% Create a figure outside the loop for all subplots
if ~isempty(subjectsForPlot) && mean(~isnan(subjectsForPlot))
    figure;
end

for subj = 1:numSubjects
    maxSNR = -inf; % Initialize maximum SNR to negative infinity
    optimalTimePointIdx = NaN; % Initialize index of optimal time point
    timePointMid = (timePoints(:,1) + timePoints(:,2)) / 2; % Midpoint of each time interval
    
    for i = 1:length(timePointMid)
        % Select the maximum SNR value between left and right at each time point
        currentMaxSNR = max(rSNR_left(subj, i), rSNR_right(subj, i));
        
        % Check if the midpoint of the current time interval is within the desired interval
        if timePointMid(i) >= lowerBound && timePointMid(i) <= upperBound
            % Find the maximum SNR within the interval
            if currentMaxSNR > maxSNR
                maxSNR = currentMaxSNR;
                optimalTimePointIdx = i;
                maxSNR_opt(subj) = currentMaxSNR;
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
%             subplot(subplotSetup(1), subplotSetup(2), find(subjectsForPlot == subj));
            
            hold on;
            % Plot individual SNR lines
            p1 = plot(timePointMid, rSNR_left(subj, :), 'g-', 'DisplayName', 'rSNR Left'); % Plot left SNR in green
            p2 = plot(timePointMid, rSNR_right(subj, :), 'r-', 'DisplayName', 'rSNR Right'); % Plot right SNR in red
            
            % Highlight the optimal time point
            p4 = scatter(timePointMid(optimalTimePointIdx), maxSNR_opt(subj), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimal Point'); % Mark optimal point as a black circle
            
            % Plot vertical lines for the lower and upper bounds
            xline(lowerBound, '--k', 'DisplayName', 'Lower Bound');
            xline(upperBound, '--k', 'DisplayName', 'Upper Bound');
            
            hold off;
            set(gca, 'color', 'none');
            title(sprintf('S%d|Max SNR:%.1f', subj, maxSNR_opt(subj)));
            box off;
            set(gcf, 'Position', [200, 400, 800, 400]);
        end
    end
end

if ~isnan(subjectsForPlot)
    legend('show');
    xlabel('Time Points');
    ylabel('rSNR');
end

end
