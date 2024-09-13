function plotDiscordantReactionTimes(sub_MF_pt, discordantSubs, T_patn_MEGfMRI)
    % plotDiscordantReactionTimes: Plot reaction times for discordant subjects
    %
    % Inputs:
    %   - sub_MF_pt: List of subject IDs
    %   - discordantSubs: Indices of discordant subjects
    %   - T_patn_MEGfMRI: Table containing reaction time data
    %   - save_dir: Directory to save the plot
    %
    % Outputs:
    %   - A bar plot of reaction times for discordant subjects saved as an SVG file

    % Extract discordant subjects and their corresponding reaction times
    discordant_subs = sub_MF_pt(discordantSubs);
    discordant_indices = ismember(T_patn_MEGfMRI.Sub_ID, discordant_subs);
    discordant_RT = T_patn_MEGfMRI.Avg(discordant_indices);
    meanRT = nanmean(discordant_RT);

    % Prepare x-axis labels for the bar plot
    xllabel = arrayfun(@(x) ['Subj ', num2str(x)], discordantSubs, 'UniformOutput', false);

    % Plot the reaction times of discordant subjects
    figure;
    bar(discordant_RT);
    xlabel('Subjects');
    set(gca, 'color', 'none');
    ylabel('Response Time (sec)');
    title({'Response Time', 'of Discordant Subjects'});
    set(gca, 'XTick', 1:length(discordant_subs), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
    set(gcf, 'Position', [1000, 100, 300, 300]);
    hold on;
    line(get(gca, 'xlim'), [meanRT meanRT], 'Color', 'red', 'LineStyle', '--');
    hold off;

%     % Save the figure
%     cfg = [];
%     cfg.outdir = save_dir;
%     filename = 'ReactionTime_Discordant_LIs';
%     cfg.filename = filename;
%     cfg.type = 'svg';
%     do_export_fig(cfg);
%     close all;
%     combined_path = fullfile(save_dir, [cfg.filename, '.svg']);
%     web(combined_path, '-new');
end
