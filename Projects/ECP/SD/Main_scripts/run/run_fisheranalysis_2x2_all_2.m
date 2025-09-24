
%% Summaries for multiple measures, now with a "Comparison" column
% clc, close all;

% Initialize an empty summary table
% Tsummary = table('Size',[0 5], ...
%     'VariableTypes',["string","string","string","double","double"], ...
%     'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});

%%
% Tsummary = table( ...
%     'Size',[0 11], ...
%     'VariableTypes',["string","string","string", "double","double","double","double", ...
%     "double","double","double","double"], ...
%     'VariableNames',{'Measure','Comparison','ROI', ...
%     'pVal','OR','CI_low','CI_high', ...
%     'A_disc','B_disc','A_conc','B_conc'});



% varNames = {'Measure','Comparison','ROI', ...
%             'qVal','OR','CI','obsCounts'};   % ? no pVal column

% Tsummary = table( ...
%     'Size',[0 numel(varNames)], ...
%     'VariableTypes',["string","string","string", ...
%                      "double","double","double","double"], ...
%     'VariableNames',varNames);

% Tsummary = [];

% Initialize an empty summary table
% Tsummary = table('Size',[0 5], ...
%     'VariableTypes',["string","string","string","double","double"], ...
%     'VariableNames',{'Measure','Comparison','ROI','qVal','OR'});
% 
% Tsummary = table( ...
%     'Size',[0 7], ...
%     'VariableTypes',["string","string","string", ...
%     "double","double","double","double"], ...
%     'VariableNames',{'Measure','Comparison','ROI', ...
%     'qVals','OR','CI','obsCounts'});

% --- master summary table -----------------------------------------
Tsummary = table( ...
    'Size',[0 7], ...
    'VariableTypes',["string","string","string", ...
                     "double","double","double","double"], ...
    'VariableNames',{'Measure','Comparison','ROI', ...
                     'qVal','OR','CI','obsCounts'});   % ? singular


disp('Fisheranalysis 2x2 (Chi-Square)')

alphaLevel = 0.05;                        % FDR threshold

%% 1) Symbol ACC
cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical    = T1_epil_measures_upted.SymbolACCcat;
cfg.discordColumn    = 'Gross_Discord_Subs';
cfg.categoryList     = {'Low','Mid','High'};
cfg.doFisher         = true;
cfg.binA             = {'Low'};
cfg.binB             = {'High', 'Mid'};
cfg.title            = "Symbol ACC";
cfg.nSubjects        = length(T1_epil_measures_upted.SubjectID);
[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg);
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 2) Animal ACC
cfg.myCategorical    = T1_epil_measures_upted.AnimalACCcat;
cfg.title            = "Animal ACC";

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

% needs Statistics Toolbox
ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 3) Animal RT
cfg.myCategorical    = T1_epil_measures_upted.AnimalRTcat;
cfg.categoryList     = {'Fast','Moderate','Slow'};
cfg.binA             = {'Moderate','Slow'};
cfg.binB             = {'Fast'};
cfg.title            = "Animal RT";

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 4) Symbol RT
cfg.myCategorical    = T1_epil_measures_upted.SymbolRTcat;
cfg.title            = "Symbol RT";

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);


ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 5) EHQ
cfg.myCategorical    = T1_epil_measures_upted.EHQcat;
cfg.categoryList     = {'Left','Right','Ambi'};
% cfg.binA             = {'Ambi','Left'};
cfg.binA             = {'Left'};
cfg.binB             = {'Right'};
cfg.title            = "EHQ";

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 6) TLE side
cfg.myCategorical    = T1_epil_measures_upted.TLEside;
cfg.categoryList     = {'Left','Right','Bilateral'};
% cfg.binA             = {'Bilateral','Right'};
cfg.binA             = {'Right'};
cfg.binB             = {'Left'};
cfg.title            = "TLE side";

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];

%% 7) AED
cfg.myCategorical    = T1_epil_measures_upted.AEDcat;
cfg.categoryList     = {'1','2','3plus'};
cfg.binA             = {'1','2'};
cfg.binB             = {'3plus'};
cfg.title            = "AED";

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 8) LTGTC
cfg.myCategorical    = T1_epil_measures_upted.LTGTCcat;
cfg.categoryList     = {'0','1-5','6-20','21plus'};
cfg.binA             = {'0','1-5','6-20'};
cfg.binB             = {'21plus'};
cfg.title            = "LTGTC";

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 9) SG
cfg.myCategorical    = T1_epil_measures_upted.SGcat;
cfg.categoryList     = {'1to2','0','3plus'};
cfg.binA             = {'1to2','0'};
cfg.binB             = {'3plus'};
cfg.title            = "SG";

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 10) CP freq
cfg.myCategorical    = T1_epil_measures_upted.cp_freq_cat;
cfg.categoryList     = {'0', '1to5','11plus','6to10'};
cfg.binA             = {'1to5', '0'};
cfg.binB             = {'11plus','6to10'};
cfg.title            = "CP freq";

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 11) tSSS_2cat
cfg = [];
cfg.bestResultsTable = bestResultsTable;  % same table used above
cfg.myCategorical    = T1_epil_measures_upted.tSSS_2cat_broad;
cfg.discordColumn    = 'Gross_Discord_Subs';   % your existing mismatch variable
cfg.categoryList     = {'Low','High'};         % the categories in tSSS_2cat
cfg.doFisher         = true;
cfg.binA             = {'Low'};                % "Group A"
cfg.binB             = {'High'};               % "Group B"
cfg.title            = 'tSSS-2cat - Broad';
cfg.nSubjects        = length(T1_epil_measures_upted.SubjectID);

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', ...
    strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_', cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 12) megnet_2cat
cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical    = T1_epil_measures_upted.megnet_2cat_broad;
cfg.discordColumn    = 'Gross_Discord_Subs';
cfg.categoryList     = {'Low','High'};
cfg.doFisher         = true;
cfg.binA             = {'Low'};
cfg.binB             = {'High'};
cfg.title            = 'megnet-2cat - Broad';
cfg.nSubjects        = length(T1_epil_measures_upted.SubjectID);

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', ...
    strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_', cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 11) tSSS_2cat
cfg = [];
cfg.bestResultsTable = bestResultsTable;  % same table used above
cfg.myCategorical    = T1_epil_measures_upted.tSSS_2cat_beta;
cfg.discordColumn    = 'Gross_Discord_Subs';   % your existing mismatch variable
cfg.categoryList     = {'Low','High'};         % the categories in tSSS_2cat
cfg.doFisher         = true;
cfg.binA             = {'Low'};                % "Group A"
cfg.binB             = {'High'};               % "Group B"
cfg.title            = 'tSSS-2cat - Beta';
cfg.nSubjects        = length(T1_epil_measures_upted.SubjectID);

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', ...
    strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_', cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% 12) megnet_2cat
cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical    = T1_epil_measures_upted.megnet_2cat_beta;
cfg.discordColumn    = 'Gross_Discord_Subs';
cfg.categoryList     = {'Low','High'};
cfg.doFisher         = true;
cfg.binA             = {'Low'};
cfg.binB             = {'High'};
cfg.title            = 'megnet-2cat - Beta';
cfg.nSubjects        = length(T1_epil_measures_upted.SubjectID);

[counts, cellCounts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar2(cfg); 
qVals = mafdr(pVals,'BHFDR',true);

ComparisonStr = sprintf('{%s} vs {%s}', ...
    strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_', cfg.title]), 'svg');

Tsummary = [Tsummary; ...
    makeSummaryBlock(cfg, cellCounts, qVals, ORvals, ComparisonStr)];


%% Final summary
% disp(Tsummary);
% % Optionally save to file:
% writetable(Tsummary,fullfile(save_dir, 'Fisher_2x2.csv'));
% 
% %%
% Specify which columns need rounding
numVars = ["qVal","OR","CI"];     % adapt if you add more

for v = numVars
    Tsummary.(v) = round(Tsummary.(v), 2);      % 2-digit decimals
end

disp(Tsummary);
writetable(Tsummary, 'Fisher_2x2_full.csv');     % CSV will now show 2 d.p.



