function plot_discordance_values(d_in, discordantSubs, xllabel, Ylabel)
    mean_d_in = nanmedian(d_in);

    figure('Color','w','Position',[1000,400,450,450]);

    % Subplot for bar
    subplot(1,4,[3 4]);
    hold on;
    tmp = d_in(discordantSubs);
    for s = 1:length(tmp)
        Swarm(s, tmp(s), 'Color', [255 105 105]/256, 'DS', 'Bar', ...
                     'SPL', 0, 'EW', 0);
    end
    xlabel('Subjects');
    ylabel(Ylabel);
    set(gca, 'color', 'none');
    set(gca, 'XTick', 1:length(discordantSubs), ...
             'XTickLabel', xllabel, ...
             'XTickLabelRotation', 90);

    % SHIFT AXIS UP (to fit text)
    ax = gca;
    pos = ax.Position;  
    pos(2) = pos(2) + 0.08;   
    pos(4) = pos(4) - 0.08;   
    ax.Position = pos;

    % Possibly keep a comfortable font:
    ax.XAxis.FontSize = 9;
    ax.YAxis.FontSize = 9;

    box off; hold off;

    line(get(gca, 'xlim'), [mean_d_in, mean_d_in], 'Color', 'k', 'LineStyle', '--');

    % Subplot for box
    ax(1) = subplot(1,4,1);
    Swarm(1, d_in, [0.2, 0.6, 0.8], 'DS', 'Box');
    set(ax(1), 'YLim', [min(d_in), max(d_in)], ...
               'XLim', [.5 1.5], ...
               'XTickLabel', 'Subject', ...
               'color','none');
    xlabel('Subjects');
    ylabel(Ylabel);
    hold on;
    d_in_discordantSubs = d_in(discordantSubs);
    hScatter = scatter(ones(size(d_in_discordantSubs)), d_in_discordantSubs, ...
                       20, 'filled', 'MarkerFaceColor', [255 105 105]/256, 'MarkerEdgeColor','none');
    hold off;
    legend(hScatter, 'Discordant', 'Location', 'southoutside');
    set(gca, 'XTick', []); % no numeric ticks
    box off;
end