% close all

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
    
    [~, ~, rsnr_max] = findIndividualOptimalTimePoints_interval_rSNR4(bestResultsTable.rSNR_L{i}, bestResultsTable.rSNR_R{i}, fMRI_LI, wi, [], bounds);
    plot_max_rSNR(rsnr_max, discordSubs)
%     sgtitle(sprintf('Max rSNR Discordant ROI: %s, Method: %s', roi, method));
    if plot_option == 1
        cfg = []; cfg.outdir = save_dir; filename = sprintf('Max rSNR Discordant ROI: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    end
    
    plotDiscordantReactionTimes2(sub_MF_pt, discordSubs, T_patn_MEGfMRI);
%     sgtitle(sprintf('Reaction Time Discordant ROI: %s, Method: %s', roi, method));
    if plot_option == 1
        cfg = []; cfg.outdir = save_dir; filename = sprintf('ReactionTime Discordant ROI: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    end
    
%     plotDiscordantTaskPerformance3(sub_MF_pt(discordSubs), meanAccBySubject_Animal, meanAccBySubject_Falsefont, xllabel, save_dir, method)
    plotDiscordantTaskPerformance3(sub_MF_pt(discordSubs), meanAccBySubject_Animal, meanAccBySubject_Falsefont, discordSubs, save_dir, method)
%     sgtitle(sprintf('Task performace: %s, Method: %s', roi, method));
    if plot_option == 1
        cfg = []; cfg.outdir = save_dir; filename = sprintf('Task performace ROI2: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    end
    
    plotCorrResponseAccReactionTime(accuracyResults_updt, T_patn_MEGfMRI)
    if plot_option == 1
        cfg = []; cfg.outdir = save_dir; filename = sprintf('Corr_ReactionTvsAcc_: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    end
    
    T1 = ecpfunc_read_epil_measures();
    
end
cd(save_dir)

%%
% Accessing the SubjID field from the structure T1
subjIDs = T1.SubjID;

xllabel2 = [];
for i=1:length(discordantSubs)
    xllabel2{i} = Pt_ID{discordantSubs(i)};
    [polyout,shapeID,vertexID] = intersect(xllabel2{i},subjIDs);
    idx1(i) = vertexID;
end

% Displaying the first 10 Subject IDs
disp('First 10 Subject IDs:');
disp(subjIDs(1:10));

% Accessing and displaying the AEDCount for the subjects
aedCounts = T1.AEDCount(idx1);
disp('AED Counts for all subjects:');
disp(aedCounts);

% Accessing categorical data for TLE side
tleSide = T1.TLEside(idx1);
disp('Temporal Lobe Epilepsy side for each subject:');
disp(tleSide);

% Example of accessing and displaying data from the xllabel cell array
disp('Labels from xllabel:');
disp(xllabel);

T2 = [];
T2.SubjID = T1.SubjID(idx1);
T2.AEDCount = T1.AEDCount(idx1);
T2.TLEside  = T1.TLEside(idx1);
T2.EHQ      = T1.EHQ(idx1);
T2.LTGTC    = T1.LTGTC(idx1);
T2.LTStatus = T1.LTStatus(idx1);
T2.AEDCount = T1.AEDCount(idx1);
T2.SP_freq  = T1.SP_freq(idx1);
T2.CP_freq  = T1.CP_freq(idx1);
T2.SG_freq  = T1.SG_freq(idx1);
T2.NP1WASI_FSIQ = T1.NP1WASI_FSIQ(idx1);

plot_discordance_values(T1.NP1WASI_FSIQ, idx1, xllabel2, 'NP1WASI-FSIQ')

if plot_option == 1
    cfg = []; cfg.outdir = save_dir; filename = 'NP1WASI-FSIQ'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
end

plot_discordance_values(T1.CP_freq, idx1, xllabel2, 'CP-freq');
if plot_option == 1
    cfg = []; cfg.outdir = save_dir; filename = 'CP-freq'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
end

plot_discordance_values(T1.EHQ, idx1, xllabel2, 'EHQ');
if plot_option == 1
    cfg = []; cfg.outdir = save_dir; filename = 'EHQ'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
end

plot_discordance_values(T1.EHQ, idx1, xllabel2, 'AED Count');
if plot_option == 1
    cfg = []; cfg.outdir = save_dir; filename = 'AED Count'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
end

