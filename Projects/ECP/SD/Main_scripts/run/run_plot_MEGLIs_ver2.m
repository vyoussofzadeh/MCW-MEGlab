network_sel = [1, 2, 6, 11]; % Define the networks to include in the plot
colors = distinguishable_colors(length(network_sel)); % Generate distinct colors for each selected network

colors = [
    0.4800    0.1200    0.6600
    0.9600    0.4900         0
    0.22 0.56 0.24
    .69 .71 .17];

colors = [
    0 114 189
    217 83 25
    237 177 32
    126 47 142
    ]/256;

for i = 1:length(LI_method_label)
    
    figure; % Open a new figure window
    hold on; % Keep the plot active to add more elements
    plotHandles = gobjects(length(network_sel), 1); % Initialize array for plot handles
    
    % Loop through the selected networks
    for net_idx = 1:length(network_sel)
        current_network = network_sel(net_idx); % Current network index
        
        % Prepare data to plot
%         LI_values = squeeze(nanmean(LI_pt_val_new.(LI_method_label{i})(current_network, :, :), 2)); % Extract and average LI values for the current method and network
        LI_values = squeeze((LI_pt_val_new.(LI_method_label{i})(current_network, :, :))); % Extract and average LI values for the current method and network
        
        % Calculate mean across methods
        meanLI = nanmean(LI_values, 1);
        % Plot the averaged LI values for the current network
        plotHandles(net_idx) = plot(wi(:,1)', meanLI, 'LineWidth', 2, 'Color', colors(net_idx,:));
%         plotHandles(net_idx) = AlphaLine(wi(:,1)',LI_values, colors(net_idx,:), 'LineWidth', 1.5);
        
        % Find the maximum LI value and its corresponding time
        [maxLI, idx] = max(meanLI);
        maxTime = mean(wi(idx,:));  % Average time at the maximum LI point
        
        % Annotate the maximum value on the plot
        %     text(maxTime, maxLI, sprintf('Mx:%.2f %.2fs', maxLI, maxTime), ...
        %         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
        text(maxTime, maxLI, sprintf('%.2fs', maxTime), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
        % Draw a vertical line at the max time
        line([maxTime maxTime], ylim, 'Color', colors(net_idx,:), 'LineWidth', 1.5, 'LineStyle', '--');
        
        legendEntries{net_idx} = net_sel_mutiple_label{network_sel(net_idx)}; % Store legend entry
    end
    
    xlabel('Time (s)');
    ylabel('Laterality Index');
    title({'MEG Laterality Index'; LI_method_label{i}});
    set(gca, 'color', 'none');
    
    % Use the plot handles for the legend to ensure continuity
    legend(plotHandles,legendEntries, 'Location', 'southoutside', 'NumColumns', 4, 'Orientation', 'horizontal');
    
    box off;
    set(gcf, 'Position', [800, 400, 500, 400]);
    hold off; % Release the plot hold
    
    % Set up configuration for exporting the figure
    cfg = [];
    cfg.outdir = save_dir; % Ensure save_dir is defined and points to a valid directory path
    cfg.filename = ['MEG_Laterality_Index_' LI_method_label{i}]; % Filename without the extension
    cfg.type = figtype; % Specify the type as 'fig'
    do_export_fig(cfg); % Call the export function
    
    cd(save_dir); % Change back to the save directory
    
%     close all, combined_path = fullfile(save_dir,[cfg.filename, ['.',figtype]]); web(combined_path, '-new');
    
end
