function plotTimeIntervalData(wi, d_in, net_sel_mutiple_label, cfg_main)
% Define custom colors
customColors = [
    0.96 0.49 0
    0.22 0.56 0.24
    0.69 0.71 0.17
    0.48 0.12 0.66
    ];

% customColors = [
%     31, 78, 121;
%     132, 60, 12;
%     127, 96, 0;
%     112, 48, 160;
%     84, 130, 53]/ 255;

% Calculate midpoints for x-axis
midpoints = mean(wi, 2);

% Initialize and configure figure
figure;
hold on;
set(gca, 'color', 'none');

% Plot data using midpoints
for i = 1:size(d_in, 1)
    plot(midpoints, d_in(i,:), 'LineWidth', 3, 'Color', customColors(i,:));
end

% Set x-ticks and labels for every other interval
xticks(midpoints(1:2:end));  % x-ticks for every other midpoint
xLabels = arrayfun(@(i) sprintf('%.1f - %.1f', wi(i, 1), wi(i, 2)), 1:length(midpoints), 'UniformOutput', false);
xticklabels(xLabels(1:2:end));  % x-labels for every other interval
xtickangle(30);  % Rotate labels by 45 degrees

% Set font size for x-axis labels
ax = gca; % Get current axes
ax.XAxis.FontSize = 8; % Set font size for x-axis only
ax.YAxis.FontSize = 8; % Set font size for y-axis only

% Additional plot settings
ylabel('LIs corr (MEG vs. fMRI)');
xlabel('Time (sec)');
if isfield(cfg_main, 'title')
    title(cfg_main.title);
end
if isfield(cfg_main, 'legend') && ~isempty(cfg_main.legend)
    legend(cfg_main.legend (cfg_main.net_sel_id), 'Location', 'southoutside', 'NumColumns', 5);
else
    legend(net_sel_mutiple_label(cfg_main.net_sel_id), 'Location', 'southoutside', 'NumColumns', 5);
end
box off;
axis tight;
set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.05));  % Adjust to prevent clipping

hold off;
end
