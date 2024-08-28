function plotDiscordantReactionTimes2(sub_MF_pt, discordantSubs, T_patn_MEGfMRI)
% plotDiscordantReactionTimes2: Plot reaction times for discordant subjects
%
% Inputs:
%   - sub_MF_pt: List of subject IDs
%   - discordantSubs: Indices of discordant subjects
%   - T_patn_MEGfMRI: Table containing reaction time data
%
% Outputs:
%   - Two plots: one box plot for all subjects with discordant subjects highlighted
%     and one bar plot for the discordant subjects only

% Extract discordant subjects and their corresponding reaction times
discordant_subs = sub_MF_pt(discordantSubs);
discordant_indices = ismember(T_patn_MEGfMRI.Sub_ID, discordant_subs);
discordant_RT = T_patn_MEGfMRI.Avg(discordant_indices);
meanRT = nanmean(discordant_RT);

% Prepare x-axis labels for the bar plot
xllabel = arrayfun(@(x) ['S', num2str(x)], discordantSubs, 'UniformOutput', false);

% Plot the reaction times of discordant subjects on a box plot
figure;
subplot 121
daboxplot(T_patn_MEGfMRI.Avg, 'groups', ones(1, numel(T_patn_MEGfMRI.Avg)), 'outsymbol', 'kx', 'xtlabels', 'All Subjects', 'fill', 1);
xlabel('All Subjects');
ylabel('Response Time (sec)');
set(gca, 'FontSize', 10);
hold on;
hScatter = scatter(ones(size(discordant_RT)), discordant_RT, 'r', 'filled');
hold off;
set(gca, 'color', 'none');
% title({'Response Time'});
l = legend(hScatter, 'Discordant', 'Location', 'southout');
% legendPos = l.Position;
% legendPos(1) = legendPos(1) + 0.02;
% legendPos(2) = legendPos(2) - 0.10; % Move legend down
% l.Position = legendPos;
set(gca, 'XTick', []);
box off;
%     set(gcf, 'Position', [1000, 100, 300, 300]);

% Plot the reaction times of discordant subjects in a bar plot
%     figure;
subplot 122
bar(discordant_RT, 0.2);
xlabel('Subjects');
set(gca, 'color', 'none');
ylabel('Response Time (sec)');
% title({'Response Time', 'of Discordant Subjects'});
set(gca, 'XTick', 1:length(discordant_subs), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
set(gcf, 'Position', [1000, 100, 600, 300]);
hold on;
line(get(gca, 'xlim'), [meanRT meanRT], 'Color', 'red', 'LineStyle', '--');
hold off;
box off;

% sgtitle(sprintf('ROI: %s, Method: %s', roi, method));


end
