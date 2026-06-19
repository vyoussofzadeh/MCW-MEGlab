function plot_discordance_values_multi(d_in_cell, discordantSubs, varLabels)
% plot_discordance_values_multi
% Plots multiple covariates side by side.
% Each subplot shows all subjects as a boxplot + scatter,
% with discordant subjects highlighted.

nVars = length(d_in_cell);

if nVars ~= length(varLabels)
    error('Number of covariates must match length of varLabels.');
end

figure('Color','w','Position',[200,200,250*nVars,450]);

for i = 1:nVars

    ax = subplot(1, nVars, i);
    hold on;

    % Get data
    data_i = d_in_cell{i};
    data_i = data_i(:);

    % Remove NaNs for plotting box/scatter
    data_clean = data_i(~isnan(data_i));

    if isempty(data_clean)
        title(varLabels{i}, 'Interpreter','none');
        ylabel(varLabels{i}, 'Interpreter','none');
        text(0.5, 0.5, 'No valid data', 'HorizontalAlignment','center');
        axis off
        continue
    end

    % -----------------------------
    % 1) Box plot for all subjects
    % -----------------------------
    boxplot(data_clean, ...
        'Positions', 1, ...
        'Colors', [0.2, 0.6, 0.8], ...
        'Symbol', '');

    % Add all subject points with small jitter
    x_all = 1 + 0.06 * randn(size(data_clean));

    scatter(x_all, data_clean, 18, ...
        'filled', ...
        'MarkerFaceColor', [0.2, 0.6, 0.8], ...
        'MarkerFaceAlpha', 0.45, ...
        'MarkerEdgeColor', 'none');

    % -----------------------------
    % 2) Highlight discordant subjects
    % -----------------------------
    discordData = data_i(discordantSubs);
    discordData = discordData(~isnan(discordData));

    x_discord = 1 + 0.10 * randn(size(discordData));

    hScatter = scatter(x_discord, discordData, 30, ...
        'filled', ...
        'MarkerFaceColor', [255 105 105]/256, ...
        'MarkerEdgeColor', 'none');

    % -----------------------------
    % 3) Median line
    % -----------------------------
    medianVal = median(data_clean);

    if exist('yline', 'file') == 2 || exist('yline', 'builtin') == 5
        yline(medianVal, 'k--', sprintf('Median=%.2f', medianVal), ...
            'LabelVerticalAlignment','bottom', ...
            'LabelHorizontalAlignment','center');
    else
        line([0.5 1.5], [medianVal medianVal], ...
            'Color', 'k', ...
            'LineStyle', '--', ...
            'LineWidth', 1);
    end

    % -----------------------------
    % Cosmetics
    % -----------------------------
    yMin = min(data_clean);
    yMax = max(data_clean);

    if yMin == yMax
        yPad = max(abs(yMin)*0.1, 1);
    else
        yPad = 0.08 * (yMax - yMin);
    end

    set(ax, ...
        'YLim', [yMin-yPad, yMax+yPad], ...
        'XLim', [0.5 1.5], ...
        'XTick', [], ...
        'Color', 'none', ...
        'FontSize', 9);

    ylabel(varLabels{i}, 'FontSize', 9, 'Interpreter','none');
    title(varLabels{i}, 'FontSize', 10, 'Interpreter','none');

    box off;

    if i == nVars
        legend(hScatter, 'Discordant', 'Location','best');
    end

    hold off;

end

end


% function plot_discordance_values_multi(d_in_cell, discordantSubs, varLabels)
% % plot_discordance_values_multi  Plots multiple covariates side by side,
% %   each with a box plot for all subjects plus a scatter overlay of
% %   discordant subjects.
% %
% %   d_in_cell:   a cell array of numeric vectors, one for each covariate
% %                (e.g. {EHQ, CP_freq, FSIQ, ...}).
% %   discordantSubs: indices (row numbers) of "discordant" subjects.
% %   varLabels:   cell array of strings for y-axis labels and subplot titles.
% 
% nVars = length(d_in_cell);
% if nVars ~= length(varLabels)
%     error('Number of covariates (d_in_cell) must match length of varLabels.');
% end
% 
% % Create a figure wide enough for all subplots
% figure('Color','w','Position',[200,200,250*nVars,450]);
% 
% for i = 1:nVars
%     ax(1) = subplot(1, nVars, i);
% 
%     % Get the data for the i-th covariate
%     data_i = d_in_cell{i};
% 
%     % 1) Box plot (all subjects)
%     %         Swarm(1, data_i, [0.2, 0.6, 0.8], 'DS', 'Box');
% 
%     Swarm(1, data_i, [0.2, 0.6, 0.8], 'DS', 'Box');
%     set(ax(1), 'YLim', [min(data_i), max(data_i)], ...
%         'XLim', [.1 1.6], ...
%         'XTickLabel', 'Subject', ...
%         'color','none');
%     hold on;
% 
%     % 2) Highlight discordant subjects
%     discordData = data_i(discordantSubs);
%     hScatter = scatter( ...
%         ones(size(discordData)), ...
%         discordData, ...
%         20, ...
%         'filled', ...
%         'MarkerFaceColor', [255 105 105]/256, ...
%         'MarkerEdgeColor', 'none');
% 
%     % Optional: Horizontal line for median
%     medianVal = nanmedian(data_i);
%     yline(medianVal, 'k--', sprintf('Median=%.2f', medianVal), ...
%         'LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','center');
% 
%     hold off;
% 
%     % Cosmetics
%     box off;
%     set(gca, 'Color', 'none', ...
%         'XTick', [], ...           % We only have 1 box per subplot
%         'FontSize', 9);
%     ylabel(varLabels{i}, 'FontSize', 9, 'Interpreter','none');
%     title(varLabels{i}, 'FontSize', 10, 'Interpreter','none');
%     if i == nVars
%         legend(hScatter, 'Discordant', 'Location','best');
%     end
% end
% end
