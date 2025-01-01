function plotDiscordantReactionTimes2(sub_MF_pt, discordantSubs, T_patn_MEGfMRI)
% plotDiscordantReactionTimes2: Plot reaction times for discordant subjects
%
% Inputs:
%   - sub_MF_pt: List of all subject IDs (numeric or string)
%   - discordantSubs: Indices (1-based) referencing the discordant subjects in sub_MF_pt
%   - T_patn_MEGfMRI: Table containing reaction time data, with fields:
%         T_patn_MEGfMRI.Sub_ID
%         T_patn_MEGfMRI.Avg
%
% Outputs:
%   - Two subplots:
%       [Left]: Box plot of response times for all subjects,
%               red dots highlight the discordant subjects.
%       [Right]: "Bar" distributions (Swarm) for each discordant subject individually.

    %% 1. Identify Discordant Subjects and Reaction Times
    % Convert numeric sub_MF_pt to string if necessary
    if isnumeric(sub_MF_pt)
        sub_MF_pt_str = arrayfun(@(x) num2str(x), sub_MF_pt, 'UniformOutput', false);
    else
        sub_MF_pt_str = sub_MF_pt;  % Already string or cellstr
    end

    % Extract IDs of discordant subjects
    discordant_subs = sub_MF_pt_str(discordantSubs);

    % Identify rows in T_patn_MEGfMRI that match discordant subject IDs
    discordant_indices = ismember(T_patn_MEGfMRI.Sub_ID, discordant_subs);

    % Get reaction times for discordant subjects
    discordant_RT = T_patn_MEGfMRI.Avg(discordant_indices);

    % For reference, calculate the median (or mean if you prefer) of discordant RT
    meanRT = nanmedian(discordant_RT);

    %% 2. Prepare x-axis labels for the 'bar' subplot
    xllabel = arrayfun(@(x) sprintf('Sub%d', x), discordantSubs, 'UniformOutput', false);

    %% 3. Create Figure and Subplots
    figure('Color','w','Position',[1000,400,450,450]);

    % ========== Subplot on the RIGHT ([3 4]) for Discordant Bars ==========
    axRight = subplot(1,4,[3 4]);
    hold(axRight,'on');

    % Plot each discordant subject as a "bar" at x=1,2,3,...
    for s = 1:length(discordant_RT)
        Swarm(s, discordant_RT(s), ...
              'Color', [0.2, 0.6, 0.8], ...
              'DS', 'Bar', ...    % DistributionStyle => Bar
              'SPL', 0, ...       % No swarm points
              'EW', 0);           % No whiskers
    end

    % Add a red dashed line for median
    yline(meanRT, '--r', 'LineWidth',1.2);

    % Format axes
    xlabel('Subjects');
    ylabel('Response Time (sec)');
    set(axRight, 'XTick', 1:length(discordant_RT), ...
                 'XTickLabel', xllabel, ...
                 'XTickLabelRotation', 45, ...
                 'Box', 'off', ...
                 'Color', 'none');

    % Shift axes upward for label clearance
    posRight = axRight.Position;
    posRight(2) = posRight(2) + 0.07;  % shift up
    posRight(4) = posRight(4) - 0.07;  % reduce height
    axRight.Position = posRight;

    hold(axRight,'off');

    % ========== Subplot on the LEFT (Index 1) for Box Plot of All Subjects ==========
    axLeft = subplot(1,4,1);
    hold(axLeft,'on');

    % Plot a single "Box" distribution at x=1 for all subjects
    allRT = T_patn_MEGfMRI.Avg;  % entire data
    Swarm(1, allRT, [0.2, 0.6, 0.8], 'DS','Box');

    % Adjust y-limits, x-limits, etc.
    set(axLeft, 'YLim',[min(allRT), max(allRT)], ...
                'XLim',[0.5 1.5], ...
                'XTick',[], ...
                'XTickLabel',{}, ...
                'Box','off', ...
                'Color','none');

    xlabel(axLeft,'All Subjects');
    ylabel(axLeft,'Response Time (sec)');

    % Overlay the discordant subjects in red at x=1
    hScatter = scatter(ones(size(discordant_RT)), discordant_RT, ...
                       20, 'r','filled');

    hold(axLeft,'off');

    % Add legend for the red points
    lgd = legend(hScatter,'Discordant','Location','southoutside','FontSize',8);
    % Manually position if needed
    lgd.Units = 'normalized';
    lgd.Position = [0.13, 0.02, 0.2, 0.04];

    % Optionally shift the left subplot to ensure it doesn't overlap figure edges
    posLeft = axLeft.Position;
    posLeft(1) = posLeft(1) + 0.02;  % shift a bit right
    axLeft.Position = posLeft;

    %% 4. Done
    % Additional styling: unify font sizes, etc.
    axLeft.FontSize = 9;
    axRight.FontSize = 9;

end


% function plotDiscordantReactionTimes2(sub_MF_pt, discordantSubs, T_patn_MEGfMRI)
% % plotDiscordantReactionTimes2: Plot reaction times for discordant subjects
% %
% % Inputs:
% %   - sub_MF_pt: List of subject IDs
% %   - discordantSubs: Indices of discordant subjects
% %   - T_patn_MEGfMRI: Table containing reaction time data
% %
% % Outputs:
% %   - Two plots: one box plot for all subjects with discordant subjects highlighted
% %     and one bar plot for the discordant subjects only
% 
% % Extract discordant subjects and their corresponding reaction times
% discordant_subs = sub_MF_pt(discordantSubs);
% discordant_indices = ismember(T_patn_MEGfMRI.Sub_ID, discordant_subs);
% discordant_RT = T_patn_MEGfMRI.Avg(discordant_indices);
% 
% % meanRT = nanmean(discordant_RT);
% meanRT = nanmedian(discordant_RT);
% 
% 
% % Prepare x-axis labels for the bar plot
% xllabel = arrayfun(@(x) ['S', num2str(x)], discordantSubs, 'UniformOutput', false);
% 
% % Plot the reaction times of discordant subjects on a box plot
% figure;
% % Plot the reaction times of discordant subjects in a bar plot
% % subplot 122
% subplot(1,4, [3 4])
% % bar(discordant_RT, 0.2);
% 
% for s = 1:length(discordant_RT)
%     Swarm(s, discordant_RT(s), 'Color', [0.2, 0.6, 0.8], 'DS', 'Bar', 'SPL', 0, 'EW', 0)
% end
% 
% xlabel('Subjects');
% set(gca, 'color', 'none');
% ylabel('Response Time (sec)');
% % title({'Response Time', 'of Discordant Subjects'});
% set(gca, 'XTick', 1:length(discordant_subs), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
% % set(gcf, 'Position', [1000, 100, 600, 300]);
% hold on;
% line(get(gca, 'xlim'), [meanRT meanRT], 'Color', 'red', 'LineStyle', '--');
% hold off;
% box off;
% 
% 
% ax(1) = subplot(1,4,1);
% % daboxplot(T_patn_MEGfMRI.Avg, 'groups', ones(1, numel(T_patn_MEGfMRI.Avg)), 'outsymbol', 'kx', 'xtlabels', 'All Subjects', 'fill', 0);
% Swarm(1, T_patn_MEGfMRI.Avg, [0.2, 0.6, 0.8], 'DS', 'Box')
% set(ax(1), 'YLim', [min(T_patn_MEGfMRI.Avg), max(T_patn_MEGfMRI.Avg)], ...
%     'XLim', [.5 1.5], ...
%     'XTickLabel', 'Subject')
% 
% xlabel('All Subjects');
% ylabel('Response Time (sec)');
% % set(gca, 'FontSize', 10);
% hold on;
% hScatter = scatter(ones(size(discordant_RT)), discordant_RT, 'r', 'filled');
% hold off;
% set(gca, 'color', 'none');
% l = legend(hScatter, 'Discordant', 'Location', 'southout');
% set(gca, 'XTick', []);
% box off;
% set(gcf, 'Position', [1000   400   450   400]);
% 
% end
