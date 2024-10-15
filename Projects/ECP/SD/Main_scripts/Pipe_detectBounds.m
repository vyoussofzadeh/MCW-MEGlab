

% Calculate the absolute difference between left and right
abs_LR = abs(rSNR_roi.left - rSNR_roi.right);

% Flatten the matrix to compute percentiles or thresholds across all elements
flattened_abs_LR = abs(abs_LR(:));

% Identify lower bound (0 for absolute difference)
lower_bound = 0;

% Identify upper bound using, for example, the 95th percentile
upper_bound = prctile(flattened_abs_LR, 95);

% Optionally, visualize the data to ensure meaningful thresholding
histogram(flattened_abs_LR);
xlabel('abs(L-R)');
ylabel('Frequency');
title('Distribution of abs(L-R)');



%%
% Initialize a new matrix to store the absolute differences
abs_LR = zeros(size(rSNR_roi.left));

% Calculate the absolute difference for each ROI
for i = 1:size(rSNR_roi.left, 1)  % Iterating over each ROI
    abs_LR(i, :, :) = abs(rSNR_roi.left(i, :, :) - rSNR_roi.right(i, :, :));
end


% Initialize vectors to store the upper bounds for each ROI
upper_bounds = zeros(1, size(rSNR_roi.left, 1));

% Calculate the upper bound for each ROI using the 95th percentile
for i = 1:size(rSNR_roi.left, 1)
    flattened_abs_LR = abs_LR(i, :, :);
    upper_bounds(i) = prctile(flattened_abs_LR(:), 95);
end


% Assume you want to analyze further or threshold based on upper_bounds
for i = 1:size(rSNR_roi.left, 1)  % Iterate over each ROI
    current_abs_LR = squeeze(abs_LR(i, :, :));  % Squeeze to remove singleton dimensions
    % Apply threshold based on upper bound
    filtered_values = current_abs_LR(current_abs_LR < upper_bounds(i));
    % Further analysis can be done here with filtered_values
end

figure,
% Example histogram of absolute differences for a selected ROI
histogram(squeeze(abs_LR(11, :, :)));
xlabel('Absolute Difference');
ylabel('Frequency');
title('Distribution of Absolute Differences for ROI 1');


%%
% Assuming abs_LR contains the absolute differences computed as before
numROIs = size(rSNR_roi.left, 1);
numSubjects = size(rSNR_roi.left, 2);
numWindows = size(rSNR_roi.left, 3);

% Matrix to store the upper bound for each window in each ROI
upper_bounds = zeros(numROIs, numWindows);

% Calculate the upper bound for each window for each ROI using a percentile
for i = 1:numROIs
    for j = 1:numWindows
        data_window = squeeze(abs_LR(i, :, j));  % Extract data for all subjects at the current window
        upper_bounds(i, j) = prctile(data_window, 95);  % Compute the 95th percentile
    end
end

% Finding best windows based on a condition, e.g., upper bounds less than a set value
threshold = 100;  % Define this based on your criteria
best_windows = upper_bounds > threshold;  % Logical index where true indicates a "best" window

clc
close all
% Example: Plotting the upper bounds and highlighting best windows
figure;
for i = 11%numROIs
    plot(upper_bounds(i, :)', 'LineWidth', 2);  % Plot upper bounds for each ROI
    hold on;
    % Highlight best windows
    plot(find(best_windows(i, :)), upper_bounds(i, best_windows(i, :)), 'ro');
end
xlabel('Window');
ylabel('Upper Bound Value');
title('Upper Bounds and Best Windows by ROI');
legend('Upper Bounds', 'Best Windows');


%%

% Assuming abs_LR contains the absolute differences computed as before
numROIs = size(rSNR_roi.left, 1);
ROIs = [1,2, 6, 11];
numSubjects = size(rSNR_roi.left, 2);
numWindows = size(rSNR_roi.left, 3);

% Matrix to store the 95th percentile for each window in each ROI
percentile95 = zeros(numROIs, numWindows);

% Calculate the 95th percentile for each window for each ROI
for i = ROIs
    for j = 1:numWindows
        data_window = squeeze(abs_LR(i, :, j));  % Extract data for all subjects at the current window
        percentile95(i, j) = prctile(data_window, 90);  % Compute the 95th percentile
    end
end

% Define the threshold as the 95th percentile
thresholds = 0.3.*max(percentile95');  % This sets a unique threshold for each window in each ROI

% Finding windows where the mean or median is below the 95th percentile threshold
best_windows = zeros(numROIs, numWindows);
for i = ROIs
    for j = 1:numWindows
        % Assuming you are comparing to the mean or median of the absolute differences
        mean_abs_diff = mean(squeeze(abs_LR(i, :, j)));
        if mean_abs_diff > thresholds(i)
            best_windows(i, j) = 1;  % Mark window as "best" if the mean is less than the 95th percentile
        end
    end
end

% Example: Plotting the 95th percentiles and highlighting the best windows
best_indices_roi = [];
clc
close all
figure;
for i = ROIs
    plot(wi(:,1)',percentile95(i, :), 'LineWidth', 2);  % Plot 95th percentile for each ROI
    hold on;
%     pause,
    % Highlight best windows
    best_indices = find(best_windows(i, :) == 1);
    best_indices = [best_indices(1), best_indices(end)];
    plot(wi(best_indices,1)', percentile95(i, best_indices), 'ro');
    best_indices_roi = [best_indices; best_indices_roi];
end
xlabel('Window');
ylabel('95th Percentile Value');
title('95th Percentile and Best Windows by ROI');
legend('95th Percentile', 'Best Windows');

disp([wi(best_indices_roi(:,1),1), wi(best_indices_roi(:,2),1)])

%%

% Initialize a matrix to store information on the best windows
best_windows = zeros(numROIs, numWindows);

% Iterate over each ROI and each window
for i = 1:numROIs
    for j = 1:numWindows
        % Calculate the mean of the absolute differences for the current window
        mean_abs_diff = mean(squeeze(abs_LR(i, :, j)));

        % Check if the mean abs(L-R) for the current window is greater than the 95th percentile threshold
        if mean_abs_diff > thresholds(i, j)
            best_windows(i, j) = 1;  % Mark window as "best" if the mean is greater than the 95th percentile
        end
    end
end



%%
% Assuming abs_LR contains the absolute differences computed as before
numROIs = size(rSNR_roi.left, 1);
numWindows = size(rSNR_roi.left, 3);

% Matrix to store the average abs(diff) for each window in each ROI
mean_abs_diff = zeros(numROIs, numWindows);

% Calculate the mean absolute difference for each window for each ROI
for i = 1:numROIs
    for j = 1:numWindows
        data_window = squeeze(abs_LR(i, :, j));  % Extract data for all subjects at the current window
        mean_abs_diff(i, j) = mean(data_window);  % Compute the mean
    end
end

% Calculate lower and upper bounds using percentiles
lower_bounds = prctile(mean_abs_diff, 85, 2);  % 25th percentile across windows for each ROI
upper_bounds = prctile(mean_abs_diff, 95, 2);  % 75th percentile across windows for each ROI

% Assuming numSubjects is the number of subjects
optimal_windows = cell(numROIs, numSubjects);

for i = 1:numROIs
    for s = 1:numSubjects
        for j = 1:numWindows
            if (abs_LR(i, s, j) >= lower_bounds(i)) && (abs_LR(i, s, j) <= upper_bounds(i))
                optimal_windows{i, s} = [optimal_windows{i, s}, j];  % Collect windows that meet the criteria
            end
        end
    end
end

% Example of plotting the results for a selected ROI and subject
selectedROI = 11;
selectedSubject = 1;
figure;
stem(optimal_windows{selectedROI, selectedSubject}, ones(size(optimal_windows{selectedROI, selectedSubject})), 'filled');
xlabel('Window');
ylabel('Selected');
title(sprintf('Optimal Windows for ROI %d, Subject %d', selectedROI, selectedSubject));


