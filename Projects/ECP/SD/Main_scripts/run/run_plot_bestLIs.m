
%% Setup / Initialization
% close all;  % Uncomment if desired to close all figures initially
% Access the SubjID field from the structure T1

%% Extract ROI and Method information from bestResultsTable
% (Currently using only row i=3; adjust as needed)
roi          = bestResultsTable.ROI{i};
method       = bestResultsTable.Best_LI_Method{i};
MEG_LI       = bestResultsTable.MEG_LI{i};
fMRI_LI      = bestResultsTable.fMRI_LI{i};

%% Prepare xllabel for each subject in discordantSubs
xllabel = cell(size(discordSubs));
for idx = 1:length(discordSubs)
    xllabel{idx} = ['S', num2str(discordSubs(idx))];
end

%% Response SNR
[~, ~, rsnr_max] = findIndividualOptimalTimePoints_interval_rSNR4( ...
    bestResultsTable.rSNR_L{i}, bestResultsTable.rSNR_R{i}, ...
    fMRI_LI, wi, [], bounds);

% Plot the max rSNR discordance
plot_discordance_values(rsnr_max, discordSubs, xllabel, 'rSNR max');

% Now export using our new utility function
doPlotExport(plot_option, save_dir, ...
    sprintf('Max rSNR Discordant ROI_%s_Method_%s', roi, method), 'svg');

%% Task Reaction Time (RT)
% Convert sub_MF_pt to a cell array of strings if numeric
if isnumeric(sub_MF_pt)
    sub_MF_pt_str = arrayfun(@num2str, sub_MF_pt, 'UniformOutput', false);
else
    sub_MF_pt_str = sub_MF_pt;  % Already string or cellstr
end

% Extract IDs of discordant subjects
discordant_subs = sub_MF_pt_str(discordSubs);

% Identify rows in T_patn_MEGfMRI that match the discordant subject IDs
discordant_indices = find(ismember(T_patn_MEGfMRI.Sub_ID, discordant_subs));

% Plot reaction times for discordant subjects
plot_discordance_values(T_patn_MEGfMRI.Avg, discordant_indices, xllabel, 'Response time (sec)');

% Export
doPlotExport(plot_option, save_dir, ...
    sprintf('ReactionTime_Discordant_ROI_%s_Method_%s', roi, method), 'svg');

%% Task Accuracy
% Convert discordantSubs like 'S12' to numeric by ignoring the 'S'
discordantSubsNumeric = cellfun(@(x) str2double(x(3:end)), discordant_subs);

% Identify Discordant Indices in each table
discIdxAnim = find(ismember(meanAccBySubject_Animal.Subject, discordantSubsNumeric));
discIdxSymb = find(ismember(meanAccBySubject_Falsefont.Subject, discordantSubsNumeric));

% Compute combined measure for each subject (Animal + Falsefont) / 2
combinedAccAll = (meanAccBySubject_Animal.mean_Animal_ACC + ...
    meanAccBySubject_Falsefont.mean_Falsefont_ACC) / 2;

% Ensure that 'xllabel' corresponds to the discordant subjects in 'discIdxAnim'
% 'xllabel' should be the same length as 'discIdxAnim'
if length(xllabel) ~= length(discIdxAnim)
    error('Length of xllabel (%d) does not match the number of discordant indices (%d).', ...
          length(xllabel), length(discIdxAnim));
end

% Plot accuracy for discordant subjects using 'xllabel'
plot_discordance_values(combinedAccAll, discIdxAnim, xllabel, 'Accuracy (%)');

% Export
doPlotExport(plot_option, save_dir, ...
    sprintf('TaskPerformance_ROI_%s_Method_%s', roi, method), 'svg');

%% Correlation: Reaction Time vs. Accuracy
plotCorrResponseAccReactionTime(accuracyResults_updt, T_patn_MEGfMRI);

% Export
doPlotExport(plot_option, save_dir, ...
    sprintf('Corr_ReactionTvsAcc_ROI_%s_Method_%s', roi, method), 'svg');

%% Other Covariates / Epilepsy Metrics
subjIDs = T1_epil_measures.SubjID;
cd(save_dir);

% Build xllabel2 and find row indices in T1 for the discordant subjects
xllabel2 = cell(size(discordSubs));
idx1     = zeros(size(discordSubs));
for jj = 1:length(discordSubs)
    xllabel2{jj} = Pt_ID{discordSubs(jj)};
    [~, ~, vertexID] = intersect(xllabel2{jj}, subjIDs);
    idx1(jj) = vertexID;
end

% Update xllabel for final display
xllabel = arrayfun(@(x) sprintf('S%d', x), discordSubs, 'UniformOutput', false);

%% Create a smaller table T2 for these discordant subjects
T2 = [];
T2.SubjID        = T1_epil_measures.SubjID(idx1);
T2.AEDCount       = T1_epil_measures.AEDCount(idx1);
T2.TLEside        = T1_epil_measures.TLEside(idx1);
T2.EHQ            = T1_epil_measures.EHQ(idx1);
T2.LTGTC          = T1_epil_measures.LTGTC(idx1);
T2.LTStatus       = T1_epil_measures.LTStatus(idx1);
T2.SP_freq        = T1_epil_measures.SP_freq(idx1);
T2.CP_freq        = T1_epil_measures.CP_freq(idx1);
T2.SG_freq        = T1_epil_measures.SG_freq(idx1);
T2.NP1WASI_FSIQ   = T1_epil_measures.NP1WASI_FSIQ(idx1);

%% Epilepsy Covariates
% 1) IQ (e.g., NP1WASI-FSIQ)
plot_discordance_values(T1_epil_measures.NP1WASI_FSIQ, idx1, xllabel, 'NP1WASI-FSIQ');
doPlotExport(plot_option, save_dir, 'NP1WASI-FSIQ', 'svg');

% 2) CP_freq
plot_discordance_values(T1_epil_measures.CP_freq, idx1, xllabel, 'CP-freq');
doPlotExport(plot_option, save_dir, 'CP-freq', 'svg');

% 3) EHQ
plot_discordance_values(T1_epil_measures.EHQ, idx1, xllabel, 'EHQ');
doPlotExport(plot_option, save_dir, 'EHQ', 'svg');

% 4) AED Count
plot_discordance_values(T1_epil_measures.AEDCount, idx1, xllabel, 'AED Count');
doPlotExport(plot_option, save_dir, 'AED_Count', 'svg');

%%



% % close all
% 
% xllabel = [];
% for i=1:length(discordantSubs)
%     xllabel{i} = ['S', num2str(discordantSubs(i))];
% end
% % Accessing the SubjID field from the structure T1
% subjIDs = T1.SubjID;
% 
% 
% %%
% % Loop through each entry in the table and plot
% % for i = 3:3%length(bestResultsTable.ROI)
% i=3;
% roi = bestResultsTable.ROI{i};
% method = bestResultsTable.Best_LI_Method{i};
% MEG_LI = bestResultsTable.MEG_LI{i};
% fMRI_LI = bestResultsTable.fMRI_LI{i};
% discordSubs = bestResultsTable.Best_Discord_Subs{i};
% 
% %% response SNR
% [~, ~, rsnr_max] = findIndividualOptimalTimePoints_interval_rSNR4(bestResultsTable.rSNR_L{i}, bestResultsTable.rSNR_R{i}, fMRI_LI, wi, [], bounds);
% %     plot_max_rSNR(rsnr_max, discordSubs)
% plot_discordance_values(rsnr_max, discordSubs, xllabel, 'rSNR max')
% 
% %     sgtitle(sprintf('Max rSNR Discordant ROI: %s, Method: %s', roi, method));
% if plot_option == 1
%     cfg = []; cfg.outdir = save_dir; filename = sprintf('Max rSNR Discordant ROI: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
% end
% 
% 
% %% Task reaction time
% if isnumeric(sub_MF_pt)
%     sub_MF_pt_str = arrayfun(@(x) num2str(x), sub_MF_pt, 'UniformOutput', false);
% else
%     sub_MF_pt_str = sub_MF_pt;  % Already string or cellstr
% end
% 
% % Extract IDs of discordant subjects
% discordant_subs = sub_MF_pt_str(discordSubs);
% % Identify rows in T_patn_MEGfMRI that match discordant subject IDs
% discordant_indices = find(ismember(T_patn_MEGfMRI.Sub_ID, discordant_subs)==1);
% 
% % Get reaction times for discordant subjects
% discordant_RT = T_patn_MEGfMRI.Avg(discordant_indices);
% 
% %     plotDiscordantReactionTimes2(sub_MF_pt, discordSubs, T_patn_MEGfMRI);
% plot_discordance_values(T_patn_MEGfMRI.Avg, discordant_indices, xllabel, 'Response time (sec)')
% 
% %     sgtitle(sprintf('Reaction Time Discordant ROI: %s, Method: %s', roi, method));
% if plot_option == 1
%     cfg = []; cfg.outdir = save_dir; filename = sprintf('ReactionTime Discordant ROI: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
% end
% 
% 
% %% Task Accuracy
% 
% % Convert Discordant Subject IDs to Numeric
% discordantSubsNumeric = cellfun(@(x) str2double(x(3:end)), discordant_subs);
% 
% % Identify Discordant Indices in Each Table
% discIdxAnim = find(ismember(meanAccBySubject_Animal.Subject, discordantSubsNumeric)==1);
% discIdxSymb = find(ismember(meanAccBySubject_Falsefont.Subject, discordantSubsNumeric)==1);
% 
% % Combined measure for each subject (Animal + Falsefont) / 2
% combinedAccAll = (meanAccBySubject_Animal.mean_Animal_ACC + ...
%     meanAccBySubject_Falsefont.mean_Falsefont_ACC) / 2;
% 
% % Discordant Subset for Combined Measure
% combinedAccDiscordant = combinedAccAll(discIdxAnim);
% 
% plot_discordance_values(combinedAccAll, discordant_indices, xllabel, 'Accuracy (%)')
% 
% % plotDiscordantTaskPerformance3(sub_MF_pt(discordSubs), meanAccBySubject_Animal, meanAccBySubject_Falsefont, discordSubs)
% %     sgtitle(sprintf('Task performace: %s, Method: %s', roi, method));
% if plot_option == 1
%     cfg = []; cfg.outdir = save_dir; filename = sprintf('Task performace ROI2: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
% end
% 
% %%
% plotCorrResponseAccReactionTime(accuracyResults_updt, T_patn_MEGfMRI)
% if plot_option == 1
%     cfg = []; cfg.outdir = save_dir; filename = sprintf('Corr_ReactionTvsAcc_: %s, Method: %s', roi, method); cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
% end
% 
% T1 = ecpfunc_read_epil_measures();
% 
% % end
% cd(save_dir)
% 
% %% Other co-variates, epilepsy metrics  
% 
% xllabel2 = [];
% idx1 = [];
% for i=1:length(discordSubs)
%     xllabel2{i} = Pt_ID{discordSubs(i)};
%     [polyout,shapeID,vertexID] = intersect(xllabel2{i},subjIDs);
%     idx1(i) = vertexID;
% end
% 
% xllabel = arrayfun(@(x) sprintf('Sub%d', x), discordSubs, 'UniformOutput', false);
% 
% % Displaying the first 10 Subject IDs
% disp('First 10 Subject IDs:');
% % disp(subjIDs(1:10));
% 
% % Accessing and displaying the AEDCount for the subjects
% aedCounts = T1.AEDCount(idx1);
% disp('AED Counts for all subjects:');
% disp(aedCounts);
% 
% % Accessing categorical data for TLE side
% tleSide = T1.TLEside(idx1);
% disp('Temporal Lobe Epilepsy side for each subject:');
% disp(tleSide);
% 
% % Example of accessing and displaying data from the xllabel cell array
% disp('Labels from xllabel:');
% disp(xllabel);
% 
% T2 = [];
% T2.SubjID = T1.SubjID(idx1);
% T2.AEDCount = T1.AEDCount(idx1);
% T2.TLEside  = T1.TLEside(idx1);
% T2.EHQ      = T1.EHQ(idx1);
% T2.LTGTC    = T1.LTGTC(idx1);
% T2.LTStatus = T1.LTStatus(idx1);
% T2.AEDCount = T1.AEDCount(idx1);
% T2.SP_freq  = T1.SP_freq(idx1);
% T2.CP_freq  = T1.CP_freq(idx1);
% T2.SG_freq  = T1.SG_freq(idx1);
% T2.NP1WASI_FSIQ = T1.NP1WASI_FSIQ(idx1);
% 
% plot_discordance_values(T1.NP1WASI_FSIQ, idx1, xllabel, 'NP1WASI-FSIQ')
% 
% if plot_option == 1
%     cfg = []; cfg.outdir = save_dir; filename = 'NP1WASI-FSIQ'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
% end
% 
% plot_discordance_values(T1.CP_freq, idx1, xllabel, 'CP-freq');
% if plot_option == 1
%     cfg = []; cfg.outdir = save_dir; filename = 'CP-freq'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
% end
% 
% plot_discordance_values(T1.EHQ, idx1, xllabel, 'EHQ');
% if plot_option == 1
%     cfg = []; cfg.outdir = save_dir; filename = 'EHQ'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
% end
% 
% plot_discordance_values(T1.AEDCount, idx1, xllabel, 'AED Count');
% if plot_option == 1
%     cfg = []; cfg.outdir = save_dir; filename = 'AED Count'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
% end
% % 
