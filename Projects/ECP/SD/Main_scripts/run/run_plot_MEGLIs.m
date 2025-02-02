network_sel = [1, 2, 6, 11]; % Define the networks to include in the plot

customColors = distinguishable_colors(length(network_sel)); % Generate distinct colors for each selected network

customColors = [
    0.96 0.49 0
    0.22 0.56 0.24
    0.69 0.71 0.17
    0.48 0.12 0.66
    ];

% customColors = [
%     0.4800    0.1200    0.6600
%     0.9600    0.4900         0
%     0.22 0.56 0.24
%     .69 .71 .17];

customColors = [
    0 114 189
    217 83 25
    237 177 32
    126 47 142
    ]/256;

figure; % Open a new figure window
hold on; % Keep the plot active to add more elements


% fig = figure('Units', 'inches', 'Position', [5, 5, 4.5, 5.5]);
% 
% %%% Amplitude discrim psychometric
% ax(1) = axes('Position', [.1 .75 .25 .225]); hold on
plotHandles = gobjects(length(network_sel), 1); % Initialize array for plot handles

% Loop through the selected networks
for net_idx = 1:length(network_sel)
    current_network = network_sel(net_idx); % Current network index

    % Prepare data to plot
    LI_values = []; % Initialize LI_values array
    for i = 1:length(LI_method_label)
        LI_values(i, :) = squeeze(nanmean(LI_pt_val_new.(LI_method_label{i})(current_network, :, :), 2)); % Extract and average LI values for the current method and network
    end

    % Calculate mean across methods
    meanLI = nanmean(LI_values, 1);
    % Plot the averaged LI values for the current network
%     plotHandles(net_idx) = plot(wi(:,1)', meanLI, 'LineWidth', 2, 'Color', colors(net_idx,:));
    plotHandles(net_idx) = AlphaLine(wi(:,1)',LI_values, customColors(net_idx,:), 'LineWidth', 1.5);
    
    % Find the maximum LI value and its corresponding time
    [maxLI, idx] = max(meanLI);
%     maxTime = mean(wi(idx,:));  % Average time at the maximum LI point
    maxTime = (wi(idx,1));  % Average time at the maximum LI point
    
    % Annotate the maximum value on the plot
    %     text(maxTime, maxLI, sprintf('Mx:%.2f %.2fs', maxLI, maxTime), ...
    %         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    text(maxTime, maxLI, sprintf('%.2fs', maxTime), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    % Draw a vertical line at the max time
    line([maxTime maxTime], ylim, 'Color', customColors(net_idx,:), 'LineWidth', 1.5, 'LineStyle', '--');
    
    legendEntries{net_idx} = net_sel_mutiple_label{network_sel(net_idx)}; % Store legend entry
end

xlabel('Time (s)');
ylabel('Laterality Index');
title({'MEG Laterality Index'; 'Over Time for Selected Networks'});
set(gca, 'color', 'none');

% Use the plot handles for the legend to ensure continuity
legend(plotHandles,legendEntries, 'Location', 'southoutside', 'NumColumns', 4, 'Orientation', 'horizontal');

box off;
set(gcf, 'Position', [800, 400, 500, 400]);
hold off; % Release the plot hold

% Set up configuration for exporting the figure
cfg = [];
cfg.outdir = save_dir; % Ensure save_dir is defined and points to a valid directory path
cfg.filename = 'MEG_Laterality_Index_Selected_Networks_Over_Time'; % Filename without the extension
cfg.type = figtype; % Specify the type as 'fig'
do_export_fig(cfg); % Call the export function

cd(save_dir); % Change back to the save directory

close all, combined_path = fullfile(save_dir,[cfg.filename, ['.',figtype]]); web(combined_path, '-new');