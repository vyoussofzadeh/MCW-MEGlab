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

%     findIndividualOptimalTimePoints_interval(MEG_LI, fMRI_LI, wi, discordSubs, lowerBound, upperBound); % NaN means no plot
%     sgtitle(sprintf('ROI: %s, Method: %s', roi, method));  
%     cfg = []; cfg.outdir = save_dir; filename = sprintf('Dicordance ROI: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    
%     findIndividualOptimalTimePoints_interval(MEG_LI, fMRI_LI, wi, 1:72, lowerBound, upperBound); % NaN means no plot
%     sgtitle(sprintf('ROI: %s, Method: %s', roi, method));
%     set(gcf, 'Position', [100, 100, 1600, 1300]);    
%     cfg = []; cfg.outdir = save_dir; filename = sprintf('indvLIs ROI: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

% rSNR_MEG = []; rSNR_MEG.rSNR_left = rSNR_left; rSNR_MEG.rSNR_right = rSNR_right;
% plotOptimalTimePointsOnMEG3(rSNR_MEG, MEG_LI, fMRI_LI, timePoints,  ...
%     mean(optimalInterval,2), discordantSubs, MEG_thre, lowerBound, upperBound);
% 
% rSNR_MEG = []; rSNR_MEG.rSNR_left = rSNR_left; rSNR_MEG.rSNR_right = rSNR_right;
% plotOptimalTimePointsOnMEG3_selective(rSNR_MEG, MEG_LI, fMRI_LI, timePoints,  ...
%     mean(optimalInterval,2), discordantSubs, MEG_thre, lowerBound, upperBound);

    plotDiscordantReactionTimes2(sub_MF_pt, discordSubs, T_patn_MEGfMRI);
    sgtitle(sprintf('Reaction Time Discordant ROI: %s, Method: %s', roi, method));
    cfg = []; cfg.outdir = save_dir; filename = sprintf('ReactionTime Discordant ROI: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    
    plotDiscordantTaskPerformance(sub_MF_pt(discordSubs), meanAccBySubject_Animal, meanAccBySubject_Falsefont, xllabel, save_dir, method)
    sgtitle(sprintf('Task performace: %s, Method: %s', roi, method));
    cfg = []; cfg.outdir = save_dir; filename = sprintf('Task performace ROI: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    
    plotDiscordantTaskPerformance2(sub_MF_pt(discordSubs), meanAccBySubject_Animal, meanAccBySubject_Falsefont, xllabel, save_dir, method)
    sgtitle(sprintf('Task performace: %s, Method: %s', roi, method));
    cfg = []; cfg.outdir = save_dir; filename = sprintf('Task performace ROI2: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    
    plotCorrResponseAccReactionTime(accuracyResults_updt, T_patn_MEGfMRI)
    cfg = []; cfg.outdir = save_dir; filename = sprintf('Corr_ReactionTvsAcc_: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

end
cd(save_dir)