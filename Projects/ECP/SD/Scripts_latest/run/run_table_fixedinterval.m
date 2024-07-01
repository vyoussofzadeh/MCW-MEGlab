% Initialize a table to store the summary
summaryTable = table();

% Loop through LI methods, metric types, and ROIs
for methodIdx = 1:length(LI_method_labels)
    for metricIdx = 1:length(metricNames)
        for roiIdx = 1:length(roi_labels)
            currentMetrics = resultsTable.Metrics(methodIdx).(metricNames{metricIdx})(roiIdx, :);
            [maxValue, maxIndex] = max(currentMetrics); % Find max value and its index
            
            % Calculate corresponding time from 'wi' using 'maxIndex'
%             maxTime = mean(wi(maxIndex, :), 2); % Assuming 'wi' contains start and end times of intervals
            Timeint = wi(maxIndex, :); % Assuming 'wi' contains start and end times of intervals
            
            % Add to the summary table
            newRow = {LI_method_labels{methodIdx}, metricNames{metricIdx}, roi_labels{roiIdx}, maxValue, Timeint};
            summaryTable = [summaryTable; newRow];
        end
    end
end

% Set column names
summaryTable.Properties.VariableNames = {'LI_Method', 'Metric_Type', 'ROI', 'Max_Value', 'Time_Interval'};

writetable(summaryTable, 'LI_Metrics_Summary.csv');

%
fid = fopen('LI_Metrics_Summary.txt', 'wt'); % Open file for writing
% Print a header
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', summaryTable.Properties.VariableNames{:});

% Loop through each row of the table and print
for i = 1:height(summaryTable)
    fprintf(fid, '%s\t%s\t%s\t%f\t%f\n', summaryTable.LI_Method{i}, summaryTable.Metric_Type{i}, ...
        summaryTable.ROI{i}, summaryTable.Max_Value(i), summaryTable.Time_Interval(i));
end

fclose(fid); % Close the file
disp('constant interval')
disp(summaryTable)