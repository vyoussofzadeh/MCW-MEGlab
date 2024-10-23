
network_sel = [1, 2, 6, 11]; % Define the networks to include in the plot
colors = distinguishable_colors(length(network_sel)); % Generate distinct colors for each selected network

for i = 1:length(LI_method_label)
    figure; % Open a new figure window
    hold on; % Keep the plot active to add more elements
    plotHandlesLeft = gobjects(length(network_sel), 1); % Initialize array for left plot handles
    plotHandlesRight = gobjects(length(network_sel), 1); % Initialize array for right plot handles
    % Loop through the selected networks
    for net_idx = 1:length(network_sel)
        current_network = network_sel(net_idx); % Current network index
        
        
        switch snr_option
            case 'thresh'
                rSNR_roi = transformPowSubTo3DArrays(rSNR_new.(LI_method_label{i})(network_sel,:));
                rSNR_left = squeeze(rSNR_roi.left(net_idx, :,:)); 
                rSNR_right = squeeze(rSNR_roi.right(net_idx, :,:));
            case 'raw'
                if i == 1
                    rSNR_roi = transformPowSubTo3DArrays_raw(rSNR_new.(LI_method_label{i})(network_sel,:));
                    rSNR_left = squeeze(rSNR_roi.left(net_idx, :,:)); 
                    rSNR_right = squeeze(rSNR_roi.right(net_idx, :,:));
                else
                    rSNR_roi = transformPowSubTo3DArrays(rSNR_new.(LI_method_label{i})(network_sel,:));
                    rSNR_left = squeeze(rSNR_roi.left(net_idx, :,:)); 
                    rSNR_right = squeeze(rSNR_roi.right(net_idx, :,:));
                end
        end
        
        %         LI_values = []; % Initialize LI_values array
        %         rSNR_roi = transformPowSubTo3DArrays(rSNR_new.(LI_method_label{i})); % Extract and average LI values for the current method and network
        %         rSNR_left = squeeze(rSNR_roi.left(network_sel(net_idx), :,:));
        %         rSNR_right = squeeze(rSNR_roi.right(network_sel(net_idx), :,:));
        
        
        % Calculate mean SNR across sessions or repeats (adjust as necessary)
        meanSNR_left = nanmean(rSNR_left, 1);
        meanSNR_right = nanmean(rSNR_right, 1);
        
        % Plot the averaged SNR values for left and right
        plotHandlesLeft(net_idx) = plot(wi(:,1)', meanSNR_left, 'LineWidth', 2, 'Color', colors(net_idx,:));
        plotHandlesRight(net_idx) = plot(wi(:,1)', meanSNR_right, 'LineWidth', 2, 'LineStyle', '--', 'Color', colors(net_idx,:));
        
        % Annotate and handle legends separately if needed
        legendEntries{net_idx} = ['Left - ' net_sel_mutiple_label{network_sel(net_idx)}]; % Store legend entry for left
        legendEntriesRight{net_idx} = ['Right - ' net_sel_mutiple_label{network_sel(net_idx)}]; % Store legend entry for right
    end
    
    xlabel('Time (s)');
    ylabel('Signal to Noise Ratio (SNR)');
    title({'MEG rSNR Index'; LI_method_label{i}});
    set(gca, 'color', 'none');
    
    % Use the plot handles for the legend to ensure continuity
    legend([plotHandlesLeft; plotHandlesRight], [legendEntries, legendEntriesRight], 'Location', 'southoutside', 'NumColumns', 4, 'Orientation', 'horizontal');
    
    box off;
    set(gcf, 'Position', [800, 400, 500, 400]);
    hold off; % Release the plot hold
    
    % Set up configuration for exporting the figure
    cfg = [];
    cfg.outdir = save_dir; % Ensure save_dir is defined and points to a valid directory path
    cfg.filename = ['rSNR' LI_method_label{i}]; % Filename without the extension
    cfg.type = figtype; % Specify the type as 'fig'
    do_export_fig(cfg); % Call the export function
    
    cd(save_dir); % Change back to the save directory
    
%     close all; 
    combined_path = fullfile(save_dir, [cfg.filename, ['.', figtype]]); 
%     web(combined_path, '-new');
    
end