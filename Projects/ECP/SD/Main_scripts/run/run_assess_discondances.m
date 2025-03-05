
%% Setup / Initialization
% close all;  % Uncomment if desired to close all figures initially
% Access the SubjID field from the structure T1

%% Extract ROI and Method information from bestResultsTable
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
% plot_discordance_values(rsnr_max, 1:72, 1:72, 'rSNR max');

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
discordant_indices = find(ismember(T_patn_MEGfMRI.SubjectID, discordant_subs));

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
subjIDs = T1_epil_measures.SubjectID;
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
T2.SubjectID        = T1_epil_measures.SubjectID(idx1);
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

