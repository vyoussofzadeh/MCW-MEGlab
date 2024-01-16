% function optimalIntervals = findIndividualOptimalIntervals(MEG_LI, fMRI_LI, timePoints, intervalSize)
%     numSubjects = size(MEG_LI, 1);
%     numIntervals = length(timePoints) - intervalSize + 1;
%     optimalIntervals = zeros(numSubjects, 2); % Start and end of optimal interval
%
%     for subj = 1:numSubjects
%         correlations = zeros(1, numIntervals);
%
%         for i = 1:numIntervals
%             intervalMEG_LI = mean(MEG_LI(subj, i:(i + intervalSize - 1)), 2);
%
%             % Check if there is variation in the data
%             if std(intervalMEG_LI) == 0 || isnan(fMRI_LI(subj))
%                 correlations(i) = NaN;
%             else
%                 corrMatrix = corrcoef(intervalMEG_LI, fMRI_LI(subj));
%                 if size(corrMatrix, 2) > 1
%                     correlations(i) = corrMatrix(1, 2);
%                 else
%                     correlations(i) = NaN; % Assign NaN if correlation cannot be computed
%                 end
%             end
%         end
%
%         [~, maxIdx] = nanmax(correlations); % Use nanmax to ignore NaN values
%         optimalIntervals(subj, :) = [timePoints(maxIdx), timePoints(min(maxIdx + intervalSize - 1, length(timePoints)))];
%     end
% end
%
% % Usage Example:
% % optimalIntervals = findIndividualOptimalIntervals(MEG_LI_Data, fMRI_LI_Data, TimePoints, IntervalSize);



% function optimalIntervals = findIndividualOptimalIntervals(MEG_LI, fMRI_LI, timePoints, intervalSize, stepSize)
%     numSubjects = size(MEG_LI, 1);
%     numIntervals = floor((length(timePoints) - intervalSize) / stepSize) + 1;
%     optimalIntervals = zeros(numSubjects, 2); % Start and end of optimal interval
%
%     for subj = 1:numSubjects
%         correlations = zeros(1, numIntervals);
%
%         for i = 1:stepSize:(numIntervals-1)*stepSize+1
%             intervalIdx = i:(i + intervalSize - 1);
%             intervalMEG_LI = MEG_LI(subj, intervalIdx);
%
%             if std(intervalMEG_LI) == 0 || isnan(fMRI_LI(subj))
%                 correlations(i) = NaN;
%             else
%                 corrMatrix = corrcoef(intervalMEG_LI, fMRI_LI(subj));
%                 if size(corrMatrix, 2) > 1
%                     correlations(i) = corrMatrix(1, 2);
%                 else
%                     correlations(i) = NaN; % Assign NaN if correlation cannot be computed
%                 end
%             end
%         end
%
%         [~, maxIdx] = nanmax(correlations); % Use nanmax to ignore NaN values
%         optimalIntervals(subj, :) = [timePoints(maxIdx), timePoints(min(maxIdx + intervalSize - 1, length(timePoints)))];
%     end
% end
%
% % Usage Example:
% % optimalIntervals = findIndividualOptimalIntervals(MEG_LI_Data, fMRI_LI_Data, TimePoints, IntervalSize, StepSize);


% function optimalIntervals = findIndividualOptimalIntervals(MEG_LI, fMRI_LI, timePoints, intervalSize, stepSize)
%     numSubjects = size(MEG_LI, 1);
%     numIntervals = floor((length(timePoints) - intervalSize) / stepSize) + 1;
%     optimalIntervals = zeros(numSubjects, 2); % Start and end of optimal interval
%
%     for subj = 1:numSubjects
%         rSquaredValues = zeros(1, numIntervals);
%
%         for i = 1:stepSize:(numIntervals-1)*stepSize+1
%             intervalIdx = i:(i + intervalSize - 1);
%             intervalMEG_LI = MEG_LI(subj, intervalIdx);
%
%             % Linear Regression
%             X = [ones(length(intervalIdx), 1), intervalMEG_LI']; % Predictor matrix
%             y = repmat(fMRI_LI(subj), length(intervalIdx), 1); % Response variable
%             b = regress(y, X); % Regression coefficients
%
%             yPred = X * b; % Predicted values
%             SSres = sum((y - yPred).^2); % Residual sum of squares
%             SStot = sum((y - mean(y)).^2); % Total sum of squares
%             rSquaredValues(i) = 1 - (SSres/SStot); % R² value
%         end
%
%         [~, maxIdx] = max(rSquaredValues);
%         optimalIntervals(subj, :) = [timePoints(maxIdx), timePoints(min(maxIdx + intervalSize - 1, length(timePoints)))];
%     end
% end
%
% % Usage Example:
% % optimalIntervals = findIndividualOptimalIntervals(MEG_LI_Data, fMRI_LI_Data, TimePoints, IntervalSize, StepSize);


% function optimalIntervals = findIndividualOptimalIntervals(MEG_LI, fMRI_LI, timePoints, intervalSize, stepSize, subjectForPlot)
% numSubjects = size(MEG_LI, 1);
% numIntervals = floor((length(timePoints) - intervalSize) / stepSize) + 1;
% optimalIntervals = zeros(numSubjects, 2); % Start and end of optimal interval
%
% for subj = 1:numSubjects
%     rSquaredValues = zeros(1, numIntervals);
%
%     for i = 1:stepSize:(numIntervals-1)*stepSize+1
%         intervalIdx = i:(i + intervalSize - 1);
%         intervalMEG_LI = MEG_LI(subj, intervalIdx);
%         disp(num2str(timePoints(intervalIdx(end)) - timePoints(intervalIdx(1))))
%
%         % Linear Regression
%         X = [ones(length(intervalIdx), 1), intervalMEG_LI']; % Predictor matrix
%         y = repmat(fMRI_LI(subj), length(intervalIdx), 1); % Response variable
%         b = regress(y, X); % Regression coefficients
%
%         yPred = X * b; % Predicted values
%         SSres = sum((y - yPred).^2); % Residual sum of squares
%         SStot = sum((y - mean(y)).^2); % Total sum of squares
%         rSquaredValues(i) = 1 - (SSres/SStot); % R² value
%     end
%
%     [~, maxIdx] = max(rSquaredValues);
%
%     optimalIntervals(subj, :) = [timePoints(maxIdx), timePoints(min(maxIdx + intervalSize - 1, length(timePoints)))];
%
%     % Plotting for the specified subject
%     if subj == subjectForPlot
%         figure;
%         plot(rSquaredValues);
%         title(sprintf('R² Values for Subject %d', subj));
%         xlabel('Interval Index');
%         ylabel('R² Value');
%         grid on;
%     end
% end
% end
%
% % Usage Example:
% % optimalIntervals = findIndividualOptimalIntervals(MEG_LI_Data, fMRI_LI_Data, TimePoints, IntervalSize, StepSize, SubjectForPlot);


% function optimalIntervals = findIndividualOptimalIntervals(MEG_LI, fMRI_LI, timePoints, intervalSize, stepSize, subjectForPlot)
% numSubjects = size(MEG_LI, 1);
% numIntervals = floor((length(timePoints) - intervalSize) / stepSize) + 1;
% optimalIntervals = zeros(numSubjects, 2); % Start and end of optimal interval
%
% for subj = 1:numSubjects
%     rSquaredValues = zeros(1, numIntervals);
%
%     for i = 1:numIntervals
%         intervalStartIdx = ((i - 1) * stepSize) + 1;
%         intervalEndIdx = intervalStartIdx + intervalSize - 1;
%
%         % Ensure the interval does not exceed the bounds of timePoints
%         if intervalEndIdx > length(timePoints)
%             intervalEndIdx = length(timePoints);
%             intervalStartIdx = intervalEndIdx - intervalSize + 1;
%         end
%
%         intervalMEG_LI = MEG_LI(subj, intervalStartIdx:intervalEndIdx);
%
%         % Linear Regression
%         X = [ones(length(intervalMEG_LI), 1), intervalMEG_LI']; % Predictor matrix
%         y = repmat(fMRI_LI(subj), length(intervalMEG_LI), 1); % Response variable
%         b = regress(y, X); % Regression coefficients
%
%         yPred = X * b; % Predicted values
%         SSres = sum((y - yPred).^2); % Residual sum of squares
%         SStot = sum((y - mean(y)).^2); % Total sum of squares
%         rSquaredValues(i) = 1 - (SSres/SStot); % R² value
%     end
%
%     [~, maxIdx] = max(rSquaredValues);
%     optimalStartIdx = ((maxIdx - 1) * stepSize) + 1;
%     optimalEndIdx = optimalStartIdx + intervalSize - 1;
%
%     % Ensure the optimal interval does not exceed the bounds of timePoints
%     if optimalEndIdx > length(timePoints)
%         optimalEndIdx = length(timePoints);
%     end
%
%     optimalIntervals(subj, :) = [timePoints(optimalStartIdx), timePoints(optimalEndIdx)];
%
%     % Plotting for the specified subject
%     if subj == subjectForPlot
%         figure;
%         plot(rSquaredValues);
%         title(sprintf('R² Values for Subject %d', subj));
%         xlabel('Interval Index');
%         ylabel('R² Value');
%         grid on;
%     end
% end
% end


% function optimalIntervals = findIndividualOptimalIntervals(MEG_LI, fMRI_LI, timePoints, intervalSize, stepSize, subjectForPlot)
% numSubjects = size(MEG_LI, 1);
% numIntervals = floor((length(timePoints) - intervalSize) / stepSize) + 1;
% optimalIntervals = zeros(numSubjects, 2); % Start and end of optimal interval
%
% for subj = 1:numSubjects
%     minDiff = inf; % Initialize minimum difference
%     optimalIdx = 1; % Initialize optimal index
%
%     for i = 1:numIntervals
%         intervalStartIdx = ((i - 1) * stepSize) + 1;
%         intervalEndIdx = intervalStartIdx + intervalSize - 1;
%
%         % Ensure the interval does not exceed the bounds of timePoints
%         if intervalEndIdx > length(timePoints)
%             intervalEndIdx = length(timePoints);
%             intervalStartIdx = intervalEndIdx - intervalSize + 1;
%         end
%
%         intervalMEG_LI = MEG_LI(subj, intervalStartIdx:intervalEndIdx);
%         meanIntervalMEG_LI = mean(intervalMEG_LI); % Mean MEG LI for the interval
%
%         % Check the difference with fMRI LI
%         diff = abs(meanIntervalMEG_LI - fMRI_LI(subj));
%         if diff < minDiff
%             minDiff = diff;
%             optimalIdx = i;
%         end
%     end
%
%     optimalStartIdx = ((optimalIdx - 1) * stepSize) + 1;
%     optimalEndIdx = optimalStartIdx + intervalSize - 1;
%
%     % Ensure the optimal interval does not exceed the bounds of timePoints
%     if optimalEndIdx > length(timePoints)
%         optimalEndIdx = length(timePoints);
%     end
%
%     optimalIntervals(subj, :) = [timePoints(optimalStartIdx), timePoints(optimalEndIdx)];
%
%     % Optional Plotting
%     % Plotting for the specified subject
%     if subj == subjectForPlot
%         figure;
%         plot(rSquaredValues);
%         title(sprintf('R² Values for Subject %d', subj));
%         xlabel('Interval Index');
%         ylabel('R² Value');
%         grid on;
%     end
% end


function optimalIntervals = findIndividualOptimalIntervals(MEG_LI, fMRI_LI, timePoints, intervalSize, stepSize, subjectForPlot)
numSubjects = size(MEG_LI, 1);
numIntervals = floor((length(timePoints) - intervalSize) / stepSize) + 1;
optimalIntervals = zeros(numSubjects, 2); % Start and end of optimal interval

for subj = 1:numSubjects
    minDiff = inf; % Initialize minimum difference
    firstMinDiffIdx = NaN; % Initialize index of first minimum difference
    
    for i = 1:numIntervals
        intervalStartIdx = ((i - 1) * stepSize) + 1;
        intervalEndIdx = intervalStartIdx + intervalSize - 1;
        
        % Ensure the interval does not exceed the bounds of timePoints
        if intervalEndIdx > length(timePoints)
            intervalEndIdx = length(timePoints);
            intervalStartIdx = intervalEndIdx - intervalSize + 1;
        end
        
        intervalMEG_LI = MEG_LI(subj, intervalStartIdx:intervalEndIdx);
        meanIntervalMEG_LI = mean(intervalMEG_LI); % Mean MEG LI for the interval
        
        % Check the difference with fMRI LI
        diff = abs(meanIntervalMEG_LI - fMRI_LI(subj));
        if diff < minDiff
            minDiff = diff;
            firstMinDiffIdx = i;
        end
    end
    
    optimalStartIdx = ((firstMinDiffIdx - 1) * stepSize) + 1;
    optimalEndIdx = optimalStartIdx + intervalSize - 1;
    
    % Ensure the optimal interval does not exceed the bounds of timePoints
    if optimalEndIdx > length(timePoints)
        optimalEndIdx = length(timePoints);
    end
    
    optimalIntervals(subj, :) = [timePoints(optimalStartIdx), timePoints(optimalEndIdx)];
    
    % Optional Plotting
    if subj == subjectForPlot
        figure,
        intervalStartIdx = find(timePoints == optimalIntervals(subj, 1));
        intervalEndIdx = find(timePoints == optimalIntervals(subj, 2));
        meanOptimalMEG_LI = mean(MEG_LI(subj, intervalStartIdx:intervalEndIdx));
        
        plot(timePoints, MEG_LI(subj, :), 'b'); % Plot MEG LI as a blue line
        hold on;
        
        % Highlight the optimal interval
        patch([timePoints(intervalStartIdx), timePoints(intervalEndIdx), timePoints(intervalEndIdx), timePoints(intervalStartIdx)], ...
            [min(MEG_LI(subj, :)), min(MEG_LI(subj, :)), max(MEG_LI(subj, :)), max(MEG_LI(subj, :))], ...
            'red', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Adjust color and transparency as needed
        
        % Set the y-axis limits
        ylim([-100 100]);
        
        hold off;
        title(sprintf('Subject %d | fMRI LI: %.2f | Mean MEG LI: %.2f', subj, fMRI_LI(subj), meanOptimalMEG_LI));
        xlabel('Time Points');
        ylabel('MEG LI');
    end
end



