function plot_discordance_values_multi(d_in_cell, discordantSubs, varLabels)
% plot_discordance_values_multi  Plots multiple covariates side by side,
%   each with a box plot for all subjects plus a scatter overlay of
%   discordant subjects.
%
%   d_in_cell:   a cell array of numeric vectors, one for each covariate
%                (e.g. {EHQ, CP_freq, FSIQ, ...}).
%   discordantSubs: indices (row numbers) of "discordant" subjects.
%   varLabels:   cell array of strings for y-axis labels and subplot titles.

nVars = length(d_in_cell);
if nVars ~= length(varLabels)
    error('Number of covariates (d_in_cell) must match length of varLabels.');
end

% Create a figure wide enough for all subplots
figure('Color','w','Position',[200,200,250*nVars,450]);

for i = 1:nVars
    ax(1) = subplot(1, nVars, i);
    
    % Get the data for the i-th covariate
    data_i = d_in_cell{i};
    
    % 1) Box plot (all subjects)
    %         Swarm(1, data_i, [0.2, 0.6, 0.8], 'DS', 'Box');
    
    Swarm(1, data_i, [0.2, 0.6, 0.8], 'DS', 'Box');
    set(ax(1), 'YLim', [min(data_i), max(data_i)], ...
        'XLim', [.1 1.6], ...
        'XTickLabel', 'Subject', ...
        'color','none');
    hold on;
    
    % 2) Highlight discordant subjects
    discordData = data_i(discordantSubs);
    hScatter = scatter( ...
        ones(size(discordData)), ...
        discordData, ...
        20, ...
        'filled', ...
        'MarkerFaceColor', [255 105 105]/256, ...
        'MarkerEdgeColor', 'none');
    
    % Optional: Horizontal line for median
    medianVal = nanmedian(data_i);
    yline(medianVal, 'k--', sprintf('Median=%.2f', medianVal), ...
        'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','center');
    
    hold off;
    
    % Cosmetics
    box off;
    set(gca, 'Color', 'none', ...
        'XTick', [], ...           % We only have 1 box per subplot
        'FontSize', 9);
    ylabel(varLabels{i}, 'FontSize', 9, 'Interpreter','none');
    title(varLabels{i}, 'FontSize', 10, 'Interpreter','none');
    if i == nVars
        legend(hScatter, 'Discordant', 'Location','best');
    end
end
end
