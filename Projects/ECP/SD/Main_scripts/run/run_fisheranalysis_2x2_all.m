
%% Summaries for multiple measures, now with a "Comparison" column
% clc, close all;

% Initialize an empty summary table
Tsummary = table('Size',[0 5], ...
    'VariableTypes',["string","string","string","double","double"], ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});

disp('Fisheranalysis 2x2 (Chi-Square)')

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
[counts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar(cfg);
ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

% Summarize in Tloc
nROIs = height(cfg.bestResultsTable);
Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

%% 2) Animal ACC
cfg.myCategorical    = T1_epil_measures_upted.AnimalACCcat;
cfg.title            = "Animal ACC";
[counts, pVals, ORvals, ~] = ecpfunc_plotDiscordantStackedBar(cfg);
ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

nROIs = height(cfg.bestResultsTable);
Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

%% 3) Animal RT
cfg.myCategorical    = T1_epil_measures_upted.AnimalRTcat;
cfg.categoryList     = {'Fast','Moderate','Slow'};
cfg.binA             = {'Moderate','Slow'};
cfg.binB             = {'Fast'};
cfg.title            = "Animal RT";
[counts, pVals, ORvals, ~] = ecpfunc_plotDiscordantStackedBar(cfg);
ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

nROIs = height(cfg.bestResultsTable);
Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

%% 4) Symbol RT
cfg.myCategorical    = T1_epil_measures_upted.SymbolRTcat;
cfg.title            = "Symbol RT";
[counts, pVals, ORvals, ~] = ecpfunc_plotDiscordantStackedBar(cfg);
ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');


nROIs = height(cfg.bestResultsTable);
Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

%% 5) EHQ
cfg.myCategorical    = T1_epil_measures_upted.EHQcat;
cfg.categoryList     = {'Left','Right','Ambi'};
% cfg.binA             = {'Ambi','Left'};
cfg.binA             = {'Left'};
cfg.binB             = {'Right'};
cfg.title            = "EHQ";
[counts, pVals, ORvals, ~] = ecpfunc_plotDiscordantStackedBar(cfg);
ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

%% 6) TLE side
cfg.myCategorical    = T1_epil_measures_upted.TLEside;
cfg.categoryList     = {'Left','Right','Bilateral'};
% cfg.binA             = {'Bilateral','Right'};
cfg.binA             = {'Right'};
cfg.binB             = {'Left'};
cfg.title            = "TLE side";
[counts, pVals, ORvals, ~] = ecpfunc_plotDiscordantStackedBar(cfg);
ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

%% 7) AED
cfg.myCategorical    = T1_epil_measures_upted.AEDcat;
cfg.categoryList     = {'1','2','3plus'};
cfg.binA             = {'1','2'};
cfg.binB             = {'3plus'};
cfg.title            = "AED";
[counts, pVals, ORvals, ~] = ecpfunc_plotDiscordantStackedBar(cfg);
ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');


Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

%% 8) LTGTC
cfg.myCategorical    = T1_epil_measures_upted.LTGTCcat;
cfg.categoryList     = {'0','1-5','6-20','21plus'};
cfg.binA             = {'0','1-5','6-20'};
cfg.binB             = {'21plus'};
cfg.title            = "LTGTC";
[counts, pVals, ORvals, ~] = ecpfunc_plotDiscordantStackedBar(cfg);
ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

%% 9) SG
cfg.myCategorical    = T1_epil_measures_upted.SGcat;
cfg.categoryList     = {'1to2','0','3plus'};
cfg.binA             = {'1to2','0'};
cfg.binB             = {'3plus'};
cfg.title            = "SG";
[counts, pVals, ORvals, ~] = ecpfunc_plotDiscordantStackedBar(cfg);
ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

%% 10) CP freq
cfg.myCategorical    = T1_epil_measures_upted.cp_freq_cat;
cfg.categoryList     = {'0', '1to5','11plus','6to10'};
cfg.binA             = {'1to5', '0'};
cfg.binB             = {'11plus','6to10'};
cfg.title            = "CP freq";
[counts, pVals, ORvals, ~] = ecpfunc_plotDiscordantStackedBar(cfg);
ComparisonStr = sprintf('{%s} vs {%s}', strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_'+ cfg.title]), 'svg');

Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

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

[counts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar(cfg);

ComparisonStr = sprintf('{%s} vs {%s}', ...
    strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_', cfg.title]), 'svg');

% Summarize in Tsummary
nROIs = height(cfg.bestResultsTable);
Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

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

[counts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar(cfg);

ComparisonStr = sprintf('{%s} vs {%s}', ...
    strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_', cfg.title]), 'svg');

% Summarize in Tsummary
nROIs = height(cfg.bestResultsTable);
Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

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

[counts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar(cfg);

ComparisonStr = sprintf('{%s} vs {%s}', ...
    strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_', cfg.title]), 'svg');

% Summarize in Tsummary
nROIs = height(cfg.bestResultsTable);
Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];

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

[counts, pVals, ORvals, chi2P_vals] = ecpfunc_plotDiscordantStackedBar(cfg);

ComparisonStr = sprintf('{%s} vs {%s}', ...
    strjoin(cfg.binA,','), strjoin(cfg.binB,','));

doPlotExport(plot_option, save_dir, char(['2by2Fisher_', cfg.title]), 'svg');

% Summarize in Tsummary
nROIs = height(cfg.bestResultsTable);
Tloc = table( ...
    repmat(cfg.title,nROIs,1), ...
    repmat(ComparisonStr,nROIs,1), ...
    cfg.bestResultsTable.ROI, ...
    pVals, ...
    ORvals, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal','OR'});
Tsummary = [Tsummary; Tloc];


%% Final summary
disp(Tsummary);
% Optionally save to file:
writetable(Tsummary,fullfile(save_dir, 'Fisher_2x2.csv'));
