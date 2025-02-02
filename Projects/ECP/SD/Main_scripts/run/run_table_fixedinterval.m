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
            newRow = {LI_method_labels{methodIdx}, metricNames{metricIdx}, roi_labels{roiIdx}, maxValue, Timeint, maxIndex};
            summaryTable = [summaryTable; newRow];
        end
    end
end

% Set column names
summaryTable.Properties.VariableNames = {'LI_Method', 'Metric_Type', 'ROI', 'Max_Value', 'Time_Interval', 'maxIndex'};

% Round numeric data to two decimal places before saving to CSV
summaryTable.Max_Value = round(summaryTable.Max_Value, 2);

% If Time_Interval is numeric and needs rounding
% for i = 1:height(summaryTable)
%     summaryTable.Time_Interval(i,1) = round(summaryTable.Time_Interval(i,1), 2);
%     summaryTable.Time_Interval(i,2) = round(summaryTable.Time_Interval(i,2), 2);
% end

% Save to CSV
writetable(summaryTable, fullfile(save_dir,'LI_Metrics_Summary_fixed.csv'));


% writetable(summaryTable, 'LI_Metrics_Summary.csv');

%
fid = fopen('LI_Metrics_Summary_fixed.txt', 'wt'); % Open file for writing
% Print a header
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', summaryTable.Properties.VariableNames{:});

% Loop through each row of the table and print
for i = 1:height(summaryTable)
    fprintf(fid, '%s\t%s\t%s\t%f\t%f\n', summaryTable.LI_Method{i}, summaryTable.Metric_Type{i}, ...
        summaryTable.ROI{i}, summaryTable.Max_Value(i), summaryTable.Time_Interval(i));
end

fclose(fid); % Close the file


% MEG_LI_int = getMEGLIForIntervals(MEG_LI, summaryTable, wi);
% 
% for i = 1:numel(MEG_LI_int)
%     fprintf('Row %d, time interval = [%.2f %.2f], #cols matched = %d\n', ...
%         i, summaryTable.Time_Interval(i,1), summaryTable.Time_Interval(i,2), size(MEG_LI_cell{i},2));
% end
% 
% summaryTable = [summaryTable, MEG_LI_int];
% summaryTable.Properties.VariableNames{'Var6'} = 'MEG_LI_int';
disp('constant interval')
disp(summaryTable)