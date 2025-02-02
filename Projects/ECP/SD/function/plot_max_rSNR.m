function plot_max_rSNR(rsnr_max, discordantSubs)

    % Calculate median (used for horizontal red line)
    mean_rsnr_max = nanmedian(rsnr_max);

    % Create a figure with white background
    figure('Color','w', 'Position', [1000, 400, 450, 450]);

    %% ============ Subplot on the RIGHT (Columns 3 & 4) ============

    axRight = subplot(1,4,[3 4]);
    hold(axRight,'on');

    % Plot a "bar" distribution for each 'discordantSub'
    tmp = rsnr_max(discordantSubs);
    for s = 1:length(tmp)
        Swarm(s, tmp(s), ...
              'Color', [0.2, 0.6, 0.8], ...
              'DS', 'Bar', ...    % DistributionStyle = Bar
              'SPL', 0, ...       % SwarmPointLimit = 0 (no swarm points)
              'EW', 0);           % No whiskers
    end

    % Add horizontal line at median
    yline(mean_rsnr_max, '--r', 'LineWidth',1.2);
    
    xllabel = arrayfun(@(x) sprintf('Sub%d', x), discordantSubs, 'UniformOutput', false);

    % X-axis labels
    set(axRight, 'XTick', 1:length(discordantSubs), ...
                 'XTickLabel', xllabel, ...
                 'XTickLabelRotation', 45);

    % Shift axes upward to accommodate rotated labels
    axPos = axRight.Position;   % [left, bottom, width, height]
    axPos(2) = axPos(2) + 0.08; % move up
    axPos(4) = axPos(4) - 0.08; % decrease height
    axRight.Position = axPos;

    xlabel(axRight, 'Subjects');
    ylabel(axRight, 'rSNR max');
    box(axRight,'off');
    set(axRight,'Color','none');
    hold(axRight,'off');

    %% ============ Subplot on the LEFT (Column 1) ============

    axLeft = subplot(1,4,1);
    hold(axLeft,'on');

    % Plot entire distribution as a box at x=1
    Swarm(1, rsnr_max, [0.2, 0.6, 0.8], 'DS', 'Box');

    % Adjust y-limits, x-limits, and label for clarity
    set(axLeft, 'YLim',[min(rsnr_max) max(rsnr_max)], ...
                'XLim',[0.5 1.5], ...
                'XTick',[], ...       % no numeric x-tick
                'XTickLabel',{}, ...
                'Box','off', ...
                'Color','none');

    xlabel(axLeft,'All Subjects');
    ylabel(axLeft,'rSNR max');

    % Highlight discordant subjects with red dots at x=1
    rsnr_max_discordantSubs = rsnr_max(discordantSubs);
    hScatter = scatter(axLeft, ones(size(rsnr_max_discordantSubs)), ...
                              rsnr_max_discordantSubs, ...
                              20, 'r','filled');
    hold(axLeft,'off');

    % Add a legend for the red scatter points
    lgd = legend(hScatter,'Discordant','Location','southoutside','FontSize',9);
    % Manually adjust legend position if needed
    lgd.Units = 'normalized';
    lgd.Position = [0.13, 0.02, 0.2, 0.04];  % [left, bottom, width, height]

    %% ============ Final Touches ============

    % (Optionally) unify font sizes
    axLeft.FontSize  = 9;
    axRight.FontSize = 9;

    % If needed, you can further adjust axLeft's position
    axPosLeft = axLeft.Position;
    axPosLeft(1) = axPosLeft(1) + 0.02;  % shift slightly right
    axLeft.Position = axPosLeft;

end


% function plot_max_rSNR(rsnr_max, discordantSubs, xllabel)
% 
% mean_rsnr_max = nanmedian(rsnr_max);
% 
% figure;
% subplot(1,4, [3 4])
% % b = bar(rsnr_max(discordantSubs),0.2);
% tmp = rsnr_max(discordantSubs);
% for s = 1:length(tmp)
%     Swarm(s, tmp(s), 'Color', [0.2, 0.6, 0.8], 'DS', 'Bar', 'SPL', 0, 'EW', 0)
% end
% 
% % b(1).FaceColor = rgb(66, 165, 245);
% xlabel('Subjects');
% set(gca, 'color', 'none');
% ylabel('rSNR max');
% % title({'mean';'(Animals & Falsefont)'});
% set(gca, 'XTick', 1:length(discordantSubs), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
% hold on; hold off; box off
% set(gcf, 'Position', [1000, 100, 300, 300]);
% line(get(gca, 'xlim'), [mean_rsnr_max mean_rsnr_max], 'Color', 'red', 'LineStyle', '--');
% 
% ax(1) = subplot(1,4,1);
% % daboxplot(rsnr_max, 'groups', ones(1, numel(rsnr_max)), 'outsymbol', 'kx', 'xtlabels', 'test', 'fill', 0);
% Swarm(1, rsnr_max, [0.2, 0.6, 0.8], 'DS', 'Box')
% % Swarm(1, rsnr_max, "Color", rgb(66, 165, 245), "DistributionStyle", "Box")
% set(ax(1), 'YLim', [min(rsnr_max), max(rsnr_max)], ...
%     'XLim', [.5 1.5], ...
%     'XTickLabel', 'Subject')
% 
% xlabel('All Subjects');
% ylabel('rSNR max');
% % set(gca, 'FontSize', 10);
% rsnr_max_discordantSubs = rsnr_max(discordantSubs);
% hold on;
% hScatter = scatter(ones(size(rsnr_max_discordantSubs)), rsnr_max_discordantSubs, 'r', 'filled');
% hold off;
% % title({'mean';'(Animals & Falsefont)'});
% set(gca, 'color', 'none');
% legend(hScatter, 'Discordant', 'Location', 'southout');
% % legendPos = l.Position;
% % legendPos(1) = legendPos(1) + 0.02;
% % l.Position = legendPos;
% set(gca, 'XTick', []);
% box off;
% set(gcf, 'Position', [1000   400   450   450]);
% set(gca,'color','none');
% 
