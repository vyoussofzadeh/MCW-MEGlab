function plotPower(cfg)

% Unpack cfg variables
power_left = cfg.power_left; % Array of power values for the left hemisphere
power_right = cfg.power_right; % Array of power values for the right hemisphere
labels = cfg.labels; % Labels for your data series
colors = cfg.colors; % Colors for each line in the plot
wi = cfg.wi; % Window intervals

% Calculate global min and max for consistent y-axis limits
globalMin = min([min(power_left(:)), min(power_right(:))]);
globalMax = max([max(power_left(:)), max(power_right(:))]);

figure,
subplot 131
hold on;
for j=1:length(labels)
    plot(mean(wi'), power_left(j, :), 'LineWidth', 3, 'Color', colors(j, :));
end
set(gca, 'color', 'none', 'YLim', [globalMin globalMax]);
title(cfg.title(1)), ylabel('Power');

subplot 132
hold on;
for j=1:length(labels)
    plot(mean(wi'), power_right(j, :), 'LineWidth', 3, 'Color', colors(j, :));
end
set(gca, 'color', 'none', 'YLim', [globalMin globalMax]);
title(cfg.title(2)), ylabel('Power');

subplot 133
hold on;
for j=1:length(labels)
    plot(mean(wi'), power_left(j, :) - power_right(j, :), 'LineWidth', 3, 'Color', colors(j, :));
end
lgd = legend(labels);
set(lgd, 'Box', 'off');
set(gca, 'color', 'none');
title('left - right h'), ylabel('Power diff.');
legend('Location', 'northeast');
xlabel('Time (s)');
set(gcf, 'Position', [700, 900, 800, 700]);
% set(gca, 'color', 'none', 'YLim', [globalMin globalMax]);

end
