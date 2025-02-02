function ecpfunc_assess_discondances(cfg_main)

%% Setup / Initialization
roi_sel = cfg_main.roi_sel;
bestResultsTable = cfg_main.bestResultsTable;

discordSubs  = bestResultsTable.Best_Discord_Subs{roi_sel};
roi          = bestResultsTable.ROI{roi_sel};
method       = bestResultsTable.Best_LI_Method{roi_sel};
MEG_LI       = bestResultsTable.MEG_LI{roi_sel};
fMRI_LI      = bestResultsTable.fMRI_LI{roi_sel};
rSNR_L       = bestResultsTable.rSNR_L{roi_sel};
rSNR_R       = bestResultsTable.rSNR_R{roi_sel};
optMEG_LI    = bestResultsTable.optMEG_LI{roi_sel};

wi = cfg_main.wi;
bounds = cfg_main.bounds;
plot_option = cfg_main.plot_option;
save_dir = cfg_main.save_dir;
final_combined_updt = cfg_main.final_combined_updt;
T1_epil_measures = cfg_main.T1_epil_measures;
Pt_ID = final_combined_updt.SubjectID;

Animal_ACC = final_combined_updt.Animal_ACC;
Symbol_ACC = final_combined_updt.Symbol_ACC;

Animal_RT = final_combined_updt.Animal_RT;
Symbol_RT = final_combined_updt.Symbol_RT;


%% Prepare xllabel for each subject in discordantSubs
xllabel = cell(size(discordSubs));
for idx = 1:length(discordSubs)
    xllabel{idx} = ['S', num2str(discordSubs(idx))];
end

%% Response SNR
[~, ~, rsnr_max] = findIndividualOptimalTimePoints_interval_rSNR4( ...
    rSNR_L, rSNR_R, ...
    fMRI_LI, wi, [], bounds);

% Plot the max rSNR discordance
plot_discordance_values(rsnr_max, discordSubs, xllabel, 'rSNR max');
% plot_discordance_values(rsnr_max, 1:72, 1:72, 'rSNR max');

% Now export using our new utility function
doPlotExport(plot_option, save_dir, ...
    sprintf('Max rSNR Discordant ROI_%s_Method_%s', roi, method), 'svg');

%%
% plot_discordance_values(optMEG_LI - fMRI_LI , discordSubs, xllabel, 'MEG fMRI LI diff');

%% Plot LI & rSNR
subIDs    = 1:size(fMRI_LI,1);  % or however you store subject labels
plotMEGandRsnrDiscordant(optMEG_LI, rSNR_L, rSNR_R, discordSubs, 1:size(MEG_LI,1));
% Export
doPlotExport(plot_option, save_dir, ...
    sprintf('optMEG_LI_rSNR__%s_Method_%s', roi, method), 'svg');

% plotMEGLITimeseriesWithRSNR(MEG_LI, rSNR_L, rSNR_R, wi(:,1), discordSubs, subIDs);

%% Opt-MEG LIs vs. fMRI LIs
plotMEGvsfMRI(optMEG_LI, fMRI_LI, discordSubs, subIDs);

% Export
doPlotExport(plot_option, save_dir, ...
    sprintf('MEG_vs_fMRI LI_%s_Method_%s', roi, method), 'svg');

%% Task Reaction Time (RT)
combinedRTAll = (final_combined_updt.Animal_RT  + final_combined_updt.Symbol_RT)./2;

% Plot reaction times for discordant subjects
plot_discordance_values(combinedRTAll, discordSubs, xllabel, 'Response time (sec)');

% Export
doPlotExport(plot_option, save_dir, ...
    sprintf('ReactionTime_Discordant_ROI_%s_Method_%s', roi, method), 'svg');

%% Task Accuracy
combinedAccAll = (final_combined_updt.Animal_ACC  + final_combined_updt.Symbol_ACC)./2;

% Plot accuracy for discordant subjects using 'xllabel'
plot_discordance_values(combinedAccAll, discordSubs, xllabel, 'Accuracy (%)');

% Export
doPlotExport(plot_option, save_dir, ...
    sprintf('TaskPerformance_ROI_%s_Method_%s', roi, method), 'svg');

%% Correlation: Reaction Time vs. Accuracy
ACC = [];
ACC.animalAcc = Animal_ACC;
ACC.symbolAcc = Symbol_ACC;

RT = [];
RT.animalRT = Animal_RT;
RT.symbolRT = Symbol_RT;

plotCorrResponseAccReactionTime(ACC, RT, discordSubs);

% Export
doPlotExport(plot_option, save_dir, ...
    sprintf('Corr_ReactionTvsAcc_ROI_%s_Method_%s', roi, method), 'svg');

%% Epilepsy Metrics
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

%% Create a table of discordant subjects
T2 = [];
T2.SubjectID      = T1_epil_measures.SubjectID(idx1);
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

%% Fisher exact test
% In your ecpfunc_assess_discondances function, after you've defined 'discordSubs':
% All subjects indices (if you haven't defined them yet):
allSubs = 1:size(MEG_LI,1); % or however many subjects you have

% Non-discordant subset
nonDiscordSubs = setdiff(allSubs, discordSubs);

% Extract TLEside from T1_epil_measures
allTLE = T1_epil_measures.TLEside;  % Nx1 array

% Subset TLE side
discordTLE    = allTLE(discordSubs);
nonDiscordTLE = allTLE(nonDiscordSubs);

% Count 'Left' vs. 'Right' in each group
nLeft_discord    = sum(discordTLE == 'Left');
nRight_discord   = sum(discordTLE == 'Right');

nLeft_nonDiscord  = sum(nonDiscordTLE == 'Left');
nRight_nonDiscord = sum(nonDiscordTLE == 'Right');

% Build 2x2 contingency: [discord vs non-discord] x [Left vs Right]
ContTable = [nLeft_discord,   nRight_discord;
    nLeft_nonDiscord,nRight_nonDiscord];

% Run Fisher's exact test
[~, fisherP, stats] = fishertest(ContTable);

fprintf('Fisher exact test comparing TLE=Left/Right by Discordance:\n');
disp(ContTable);
fprintf('p-value = %.4g\n', fisherP);

% (Optional) if needed, see the odds ratio:
if isfield(stats,'OddsRatio')
    fprintf('Odds Ratio = %.3f\n', stats.OddsRatio);
    disp('95% CI:'); disp(stats.ConfidenceInterval);
end

%% Regression model

% Example: pick some columns from final_combined_updt
T = table( ...
    optMEG_LI, ...
    optMEG_LI - fMRI_LI, ...
    final_combined_updt.Animal_RT, ...
    final_combined_updt.Symbol_RT, ...
    final_combined_updt.Animal_ACC, ...
    final_combined_updt.Symbol_ACC, ...
    final_combined_updt.AEDCount, ...
    final_combined_updt.TLEside, ...
    final_combined_updt.EHQ, ...
    final_combined_updt.NP1WASI_FSIQ, ...
    final_combined_updt.CP_freq, ...
    final_combined_updt.LTGTC, ...
    final_combined_updt.NP1WASI_FSIQ, ...
    rsnr_max', ...
    fMRI_LI, ...
    'VariableNames', { ...
    'optMEG_LI', ...         % The outcome
    'MEG_fMRI_diff', ...
    'Animal_RT', ...
    'Symbol_RT', ...
    'Animal_ACC', ...
    'Symbol_ACC', ...
    'AEDCount', ...
    'TLEside', ...
    'EHQ', ...
    'FSIQ' ...            % renamed NP1WASI_FSIQ -> FSIQ for convenience
    'CP_freq', ...
    'LTGTC',...
    'NP1WASI_FSIQ', ...
    'rSNR', ...
    'fMRI_LI'
    });

% 3) Convert TLEside to categorical if necessary
if ~iscategorical(T.TLEside)
    T.TLEside = categorical(T.TLEside);
end

% (Optional) Remove rows with any missing values
T(any(ismissing(T),2), :) = [];

% 4) Fit the linear model
lm = fitlm(T, 'optMEG_LI ~ Animal_RT + Symbol_RT + Animal_ACC + Symbol_ACC + AEDCount + TLEside + EHQ + FSIQ + CP_freq + NP1WASI_FSIQ + rSNR + fMRI_LI');
disp(lm);

lm = fitlm(T, 'MEG_fMRI_diff ~ Animal_RT + Symbol_RT + Animal_ACC + Symbol_ACC + AEDCount + TLEside + EHQ + FSIQ + CP_freq + NP1WASI_FSIQ + rSNR + fMRI_LI');
disp(lm);

% 5) Inspect model output:
% - Look at 'lm.Coefficients' to see each predictor's slope, p-value
% - Check 'lm.Rsquared.Ordinary' or 'lm.Rsquared.Adjusted' for model fit

%%
figure('Color','w','Position',[1000,400,500,500]);

% 1) Subplot #1
subplot(2,1,1);
gscatter(T.EHQ, T.optMEG_LI, T.TLEside);
xlabel('EHQ');
ylabel('MEG LI');
legend('off');  % Hide legend in this subplot

% 2) Subplot #2
subplot(2,1,2);
gscatter(T.Animal_ACC, T.optMEG_LI, T.TLEside);
xlabel('Animal-ACC');
ylabel('MEG LI');
legend('off');  % Hide legend in this subplot

% 3) Subplot #3 (Show the legend here)
% subplot(3,1,3);
% gscatter(T.EHQ, T.fMRI_LI, T.TLEside);
% xlabel('EHQ');
% ylabel('fMRI LI');
legend('Location','best');  % Only show the legend in the last panel
doPlotExport(plot_option, save_dir, 'LI_vs_EHQ_TLEside', 'svg');

%%
d_in_cell = {
    T1_epil_measures.EHQ, 
    T1_epil_measures.CP_freq, 
    T1_epil_measures.NP1WASI_FSIQ
    T1_epil_measures.AEDCount
    };

% Create labels:
varLabels = {'EHQ','CP-freq','FSIQ','AEDCount'};

% Now call the function:
plot_discordance_values_multi(d_in_cell, discordSubs, varLabels);
doPlotExport(plot_option, save_dir, 'Epil_metrics', 'svg');


%%
% 1) Combine your data into a cell array (one entry per metric):
d_in_cell = {
    Animal_ACC, ...
    Symbol_ACC, ...
    Animal_RT,  ...
    Symbol_RT,  ...
    rsnr_max, ...
};

% 2) Provide labels for each subplot:
varLabels = {
    'Animal ACC (%)', ...
    'Symbol ACC (%)', ...
    'Animal RT (sec)', ...
    'Symbol RT (sec)', ...
    'rSNR max'
};

% 3) Call the function:
plot_discordance_values_multi(d_in_cell, discordSubs, varLabels);
doPlotExport(plot_option, save_dir, 'Taskperformace', 'svg');


%%
% run_plotmatrix
% 
% %%
% % % Example variables:
% numericVars = {'optMEG_LI','fMRI_LI','rSNR','Animal_RT','Symbol_RT',...
%     'Animal_ACC','Symbol_ACC','AEDCount','EHQ','CP_freq','FSIQ'};
% 
% numericVars_label = {'MEG-LI','fMRI-LI','rSNR','Anim-RT','Symb-RT',...
%     'Anim-ACC','Symb-ACC','AEDCount','EHQ','CP-freq','FSIQ'};
% 
% 
% numericVars = {'optMEG_LI','fMRI_LI','rSNR','Animal_ACC'};
% 
% numericVars_label = {'MEG-LI','fMRI-LI','rSNR','Anim-ACC'};
% 
% % Suppose T is your table, numericVars is a cell array of your variable names,
% % and numericVars_label are the corresponding axis labels.
% 
% % Extract numeric data
% X = T{:, numericVars};
% nVars = size(X, 2);
% 
% % 1) Compute pairwise correlations and p-values (e.g., Pearson)
% Corr = zeros(nVars, nVars);
% Pvals = ones(nVars, nVars);
% 
% for i = 1:nVars
%     for j = 1:nVars
%         if i ~= j
%             [R, P] = corr(X(:,i), X(:,j), 'Rows','complete');  % or 'Type','Spearman'
%             Corr(i,j)  = R;
%             Pvals(i,j) = P;
%         else
%             Corr(i,j)  = 1;   % diagonal (self-correlation)
%             Pvals(i,j) = 0;
%         end
%     end
% end
% 
% % 2) Create scatter matrix
% [h, AX] = plotmatrix(X);
% 
% % 3) Label axes
% for i = 1:nVars
%     % x-labels on bottom row
%     AX(nVars, i).XLabel.String = numericVars_label{i};
%     % y-labels on left column
%     AX(i, 1).YLabel.String = numericVars_label{i};
% end
% 
% % 4) Zero lines on off-diagonal subplots (optional)
% for row = 1:nVars
%     for col = 1:nVars
%         if row ~= col
%             ax = AX(row, col);
%             hold(ax, 'on');
%             xline(ax, 0,'k--','HandleVisibility','off');
%             yline(ax, 0,'k--','HandleVisibility','off');
%         end
%     end
% end
% 
% % 5) Reduce Font Size
% for row = 1:nVars
%     for col = 1:nVars
%         axHandle = AX(row, col);
%         axHandle.FontSize = 8;
%         axHandle.XLabel.FontSize = 8;
%         axHandle.YLabel.FontSize = 8;
%     end
% end
% 
% % 6) Indicate significant relationships on off-diagonal plots
% alphaLevel = 0.05;   % significance threshold
% for row = 1:nVars
%     for col = 1:nVars
%         if row ~= col
%             ax = AX(row,col);
%             hold(ax,'on');
%             
%             % Check significance
%             if Pvals(row,col) < alphaLevel
%                 % Add small text in the upper-left corner or specify a location
%                 text(ax, ...
%                     0.05*range(ax.XLim)+ax.XLim(1), ...  % X coordinate near left
%                     0.9*range(ax.YLim)+ax.YLim(1), ...   % Y coordinate near top
%                     sprintf('p=%.3g', Pvals(row,col)), ...
%                     'FontSize',8, 'Color','m', 'FontWeight','bold');
%             end
%         end
%     end
% end
% 
% % 7) Figure styling
% set(gcf,'Name','Pairwise Scatter Matrix with Significance','NumberTitle','off');

end

