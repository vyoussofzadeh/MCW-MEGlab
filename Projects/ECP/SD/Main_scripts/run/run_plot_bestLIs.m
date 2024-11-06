close all

xllabel = [];
for i=1:length(discordantSubs)
    xllabel{i} = ['S', num2str(discordantSubs(i))];
end

% Loop through each entry in the table and plot
for i = 3:3%length(bestResultsTable.ROI)
    
    roi = bestResultsTable.ROI{i};
    method = bestResultsTable.Best_LI_Method{i};
    MEG_LI = bestResultsTable.MEG_LI{i};
    fMRI_LI = bestResultsTable.fMRI_LI{i};
    discordSubs = bestResultsTable.Best_Discord_Subs{i};
    
%     findIndividualOptimalTimePoints_interval_rSNR4(bestResultsTable.rSNR_L{i}(discordantSubs,:), bestResultsTable.rSNR_R{i}(discordantSubs,:), fMRI_LI(discordantSubs), wi, 1:length(discordantSubs), bounds(discordantSubs,:));
    [~, ~, rsnr_max] = findIndividualOptimalTimePoints_interval_rSNR4(bestResultsTable.rSNR_L{i}, bestResultsTable.rSNR_R{i}, fMRI_LI, wi, [], bounds);
    plot_max_rSNR(rsnr_max, discordantSubs, xllabel)
    sgtitle(sprintf('Max rSNR Discordant ROI: %s, Method: %s', roi, method));
    cfg = []; cfg.outdir = save_dir; filename = sprintf('Max rSNR Discordant ROI: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    
        plotDiscordantReactionTimes2(sub_MF_pt, discordSubs, T_patn_MEGfMRI);
    sgtitle(sprintf('Reaction Time Discordant ROI: %s, Method: %s', roi, method));
    cfg = []; cfg.outdir = save_dir; filename = sprintf('ReactionTime Discordant ROI: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    
    plotDiscordantTaskPerformance3(sub_MF_pt(discordSubs), meanAccBySubject_Animal, meanAccBySubject_Falsefont, xllabel, save_dir, method)
    sgtitle(sprintf('Task performace: %s, Method: %s', roi, method));
    cfg = []; cfg.outdir = save_dir; filename = sprintf('Task performace ROI2: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    
    plotCorrResponseAccReactionTime(accuracyResults_updt, T_patn_MEGfMRI)
    cfg = []; cfg.outdir = save_dir; filename = sprintf('Corr_ReactionTvsAcc_: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    
end
cd(save_dir)