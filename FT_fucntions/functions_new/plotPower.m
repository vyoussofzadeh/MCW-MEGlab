function plotPower(cfg)

% data = cfg.data; % This could be your LI or any other metric you're plotting
power_left = cfg.power_left; % Array of power values for the left hemisphere
power_right = cfg.power_right; % Array of power values for the right hemisphere
labels = cfg.labels; % Labels for your data series
colors = cfg.colors; % Colors for each line in the plot
% titleText = cfg.titleText; % Title for the plot
wi = cfg.wi; % Window intervals

figure,
hold on;
for j=1:length(labels)
    plot(mean(wi'), power_left(j, :), 'LineWidth', 3, 'Color', colors(j, :));
end
legend([labels]);
set(gca, 'color', 'none');


figure,
hold on;
for j=1:length(labels)
    plot(mean(wi'), power_right(j, :), 'LineWidth', 3, 'Color', colors(j, :));
end
legend([labels]);
set(gca, 'color', 'none');


figure,
hold on;
for j=1:length(labels)
    plot(mean(wi'), power_left(j, :) - power_right(j, :), 'LineWidth', 3, 'Color', colors(j, :));
end
legend([labels]);
set(gca, 'color', 'none');

% if length(labels) > 5
%     legend([labels, {'Power Left', 'Power Right'}], 'Location', 'southoutside', 'NumColumns', 3);
% end

% title(titleText);
ylabel('Metric Value');
xlabel('Time (s)');
set(gca, 'color', 'none');

end


% function plotPower(pow_sub, wi)
% % Plots the average power for left and right hemispheres across intervals
% % pow_sub is an 11x30 struct array with fields: left, right
% % wi is the window intervals, assuming it is an Nx2 array where N is the number of intervals
%
% % Pre-allocate arrays to store average power values for plotting
% avgPowerLeft = zeros(1, size(pow_sub, 2)); % For each interval across all subjects
% avgPowerRight = zeros(1, size(pow_sub, 2));
%
% % Loop through each interval
% for i = 1:size(pow_sub, 2)
%     % Extract left and right power values for all subjects at interval i
%     leftPowers = [pow_sub(:, i).left];
%     rightPowers = [pow_sub(:, i).right];
%
%     % Compute the average power for left and right hemispheres
%     avgPowerLeft(i) = mean(leftPowers);
%     avgPowerRight(i) = mean(rightPowers);
% end
%
% % Plotting
% figure;
% hold on;
% plot(wi(:,1), avgPowerLeft, 'b-', 'LineWidth', 2); % Plot average power for left hemisphere
% plot(wi(:,1), avgPowerRight, 'r-', 'LineWidth', 2); % Plot average power for right hemisphere
% legend('Avg Power Left', 'Avg Power Right', 'Location', 'Best');
% xlabel('Time (s)');
% ylabel('Average Power');
% title('Average Power Across Subjects for Each Interval');
% grid on;
% hold off;
%
% end

