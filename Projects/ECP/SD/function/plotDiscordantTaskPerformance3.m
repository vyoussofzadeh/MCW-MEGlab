function plotDiscordantTaskPerformance3(discordant_subs, meanAccBySubject_Animal, meanAccBySubject_Falsefont, discordSubs, ~, ~)
% plotDiscordantTaskPerformance3: Plot task performance for discordant subjects
%
% Inputs:
%   - discordant_subs: Cell array of discordant subject IDs as strings (like 'EC123')
%   - meanAccBySubject_Animal: Table with fields:
%       .Subject (numeric IDs),
%       .mean_Animal_ACC
%   - meanAccBySubject_Falsefont: Table with fields:
%       .Subject (numeric IDs),
%       .mean_Falsefont_ACC
%   - xllabel: Cell array of labels for the x-axis
%   - ~: (unused inputs)
%
% This function creates:
%   - A box plot (left subplot) of combined Animal+Falsefont accuracies for all subjects,
%     with discordant subjects highlighted in red scatter.
%   - Individual bars (right subplot) for each discordant subjects combined accuracy,
%     labeled on the x-axis.

    %% 1. Convert Discordant Subject IDs to Numeric
    %   e.g., 'EC123' -> 123
    discordantSubsNumeric = cellfun(@(x) str2double(x(3:end)), discordant_subs);

    %% 2. Identify Discordant Indices in Each Table
    %   Where .Subject matches the numeric IDs
    discIdxAnim = ismember(meanAccBySubject_Animal.Subject, discordantSubsNumeric);
    discIdxSymb = ismember(meanAccBySubject_Falsefont.Subject, discordantSubsNumeric);

    %% 3. Compute Mean / Median Accuracies
    %   If you want median instead, replace `nanmean` with `nanmedian`
    meanTP_anim = nanmean(meanAccBySubject_Animal.mean_Animal_ACC);
    meanTP_symb = nanmean(meanAccBySubject_Falsefont.mean_Falsefont_ACC);

    % Combined measure for each subject (Animal + Falsefont) / 2
    combinedAccAll = (meanAccBySubject_Animal.mean_Animal_ACC + ...
                      meanAccBySubject_Falsefont.mean_Falsefont_ACC) / 2;
    % Mean of discordant subjects specifically (Animal measure) if you want:
    % If you want a red line for the combined measure instead, adjust below
    % e.g., line(..., [nanmean(combinedAccDiscordant) nanmean(combinedAccDiscordant)], ... )

    %% 4. Discordant Subset for Combined Measure
    combinedAccDiscordant = combinedAccAll(discIdxAnim); 
    % (since discIdxAnim == discIdxSymb for the same subjects)

    %% 5. Create Figure
    figure('Color','w','Position',[1000,400,450,450]);

    %---------------------- RIGHT SUBPLOT: Discordant Bars ------------------%
    axRight = subplot(1,4,[3 4]);
    hold(axRight, 'on');

    % Plot each discordant subjects combined accuracy as a bar
    for s = 1:length(combinedAccDiscordant)
        Swarm(s, combinedAccDiscordant(s), ...
              'Color', [0.2, 0.6, 0.8], ...
              'DS', 'Bar', ...   % DistributionStyle => Bar
              'SPL', 0, ...      % No swarm points
              'EW', 0);          % No whiskers
    end

    % Suppose you want a reference line for the Animal mean accuracy:
    % (If you'd prefer the mean of the combined measure, do so)
    midian_combinedAcc = nanmedian(combinedAccAll);
    yline(midian_combinedAcc, '--r', 'Acc median', 'LineWidth',1.2);
    
    
    xllabel = arrayfun(@(x) sprintf('Sub%d', x), discordSubs, 'UniformOutput', false);

    % Format x-axis
    xlabel('Subjects');
    ylabel('Accuracy (%)');
    set(axRight, 'XTick', 1:length(combinedAccDiscordant), ...
                 'XTickLabel', xllabel, ...
                 'XTickLabelRotation', 45, ...
                 'Color','none', 'Box','off');

    % Shift axes up to accommodate rotated labels
    posRight = axRight.Position; 
    posRight(2) = posRight(2) + 0.08;
    posRight(4) = posRight(4) - 0.08;
    axRight.Position = posRight;

    hold(axRight,'off');
    
    %%
   %---------------------- LEFT SUBPLOT: All Subjects Box ------------------%
    axLeft = subplot(1,4,1);
    hold(axLeft,'on');

    % Plot combined measure for all subjects as a box at x=1
    Swarm(1, combinedAccAll, [0.2, 0.6, 0.8], 'DS','Box');

    % Adjust axes
    set(axLeft, 'YLim',[min(combinedAccAll), max(combinedAccAll)], ...
                'XLim',[0.5 1.5], ...
                'XTick',[], 'XTickLabel',{}, ...
                'Color','none', 'Box','off');

    xlabel(axLeft, 'All Subjects');
    ylabel(axLeft, 'Accuracy (%)');

    % Overlay red scatter for the discordant subset
    hScatter = scatter(axLeft, ones(size(combinedAccDiscordant)), combinedAccDiscordant, ...
            20, 'r','filled');

    hold(axLeft,'off');

    % Legend for the red scatter
    lgd = legend(hScatter, 'Discordant','Location','southoutside','FontSize',8);
    lgd.Units = 'normalized';
    lgd.Position = [0.13, 0.02, 0.2, 0.04];

    % Slightly shift left subplot to avoid figure edges
    posLeft = axLeft.Position;
    posLeft(1) = posLeft(1) + 0.02;
    axLeft.Position = posLeft;

    %---------------------- Final Aesthetics ------------------------------%
    axLeft.FontSize  = 9;
    axRight.FontSize = 9;
end


% function plotDiscordantTaskPerformance3(discordant_subs, meanAccBySubject_Animal, meanAccBySubject_Falsefont, xllabel, ~, ~)
% % plotDiscordantTaskPerformance: Plot task performance for discordant subjects
% %
% % Inputs:
% %   - discordant_subs: Cell array of discordant subject IDs as strings
% %   - meanAccBySubject_Animal: Table containing accuracy data for the animal task
% %   - meanAccBySubject_Falsefont: Table containing accuracy data for the symbol task
% %   - xllabel: Cell array of labels for the x-axis
% %   - save_dir: Directory to save the plot
% %   - LI_method_label: Cell array of LI method labels
% %   - LI_method: Index of the current LI method
% 
% % Convert discordant subjects to numeric
% discordant_subs_numeric = cellfun(@(x) str2double(x(3:end)), discordant_subs);
% discordant_indices_anim = ismember(meanAccBySubject_Animal.Subject, discordant_subs_numeric);
% discordant_indices_symb = ismember(meanAccBySubject_Falsefont.Subject, discordant_subs_numeric);
% meanTP_anim = nanmean(meanAccBySubject_Animal.mean_Animal_ACC);
% meanTP_symb = nanmean(meanAccBySubject_Falsefont.mean_Falsefont_ACC);
% 
% meanAccBySubject_Anim_Falsefont = (meanAccBySubject_Animal.mean_Animal_ACC + meanAccBySubject_Falsefont.mean_Falsefont_ACC)./2;
% 
% % Plotting Animal Task Performance
% figure;
% subplot(1,4, [3 4])
% tmp = meanAccBySubject_Anim_Falsefont(discordant_indices_anim);
% 
% for s = 1:length(tmp)
%     Swarm(s, tmp(s), 'Color', [0.2, 0.6, 0.8], 'DS', 'Bar', 'SPL', 0, 'EW', 0)
% end
% 
% xlabel('Subjects');
% set(gca, 'color', 'none');
% ylabel('Accuracy (%)');
% % title({'mean';'(Animals & Falsefont)'});
% set(gca, 'XTick', 1:length(discordant_indices_anim), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
% hold on; hold off; box off
% set(gcf, 'Position', [1000, 100, 300, 300]);
% line(get(gca, 'xlim'), [meanTP_anim meanTP_anim], 'Color', 'red', 'LineStyle', '--');
% 
% 
% %%
% ax(1) = subplot(1,4,1);
% % daboxplot(T_patn_MEGfMRI.Avg, 'groups', ones(1, numel(T_patn_MEGfMRI.Avg)), 'outsymbol', 'kx', 'xtlabels', 'All Subjects', 'fill', 0);
% Swarm(1, meanAccBySubject_Anim_Falsefont, [0.2, 0.6, 0.8], 'DS', 'Box')
% set(ax(1), 'YLim', [min(meanAccBySubject_Anim_Falsefont), max(meanAccBySubject_Anim_Falsefont)], ...
%     'XLim', [.5 1.5], ...
%     'XTickLabel', 'Subject')
% 
% xlabel('All Subjects');
% ylabel('Accuracy (%)');
% set(gca, 'FontSize', 10);
% discordant_ACC_anim = meanAccBySubject_Anim_Falsefont(discordant_indices_anim);
% hold on;
% hScatter = scatter(ones(size(discordant_ACC_anim)), discordant_ACC_anim, 'r', 'filled');
% hold off;
% % title({'mean';'(Animals & Falsefont)'});
% set(gca, 'color', 'none');
% l = legend(hScatter, 'Discordant', 'Location', 'southout');
% set(gca, 'XTick', []);
% box off;
% set(gcf, 'Position', [1000   400   450   450]);
% set(gca,'color','none');
% 
% end
