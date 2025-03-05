%% run_fisheranalysis_2x3_optimized.m

% We assume you already have:
%  bestResultsTable            (with ROI + discord columns)
%  T1_epil_measures_upted.*   (various 3-level cat variables)
%  ecpfunc_plotDiscordantStackedBar_2x3  (the function that performs the 2×3 test)

% clc; close all;

%% 1) Initialize a summary table
Tsummary2x3 = table('Size',[0 4], ...
    'VariableTypes',["string","string","string","double"], ...
    'VariableNames',{'Measure','Comparison','ROI','pVal'});

% comparisonStr = "2×3(Discord vs Concord)";

%%
disp('Fisheranalysis 2x3 (Freedman-Halton analysis)')

%% 2) Define a helper function to avoid repeating code
% run2x3 = @(myCat, catList, measureName, doFH) ...
%     run_2x3_measure(bestResultsTable, myCat, catList, measureName, doFH, comparisonStr);

run2x3 = @(myCat, catList, measureName, doFH) ...
    run_2x3_measure(bestResultsTable, myCat, catList, measureName, doFH);

%% 3) Symbol ACC
myCat   = T1_epil_measures_upted.SymbolACCcat;  % 3-level cat
catList = {'Low','Mid','High'};
measure = 'Symbol ACC (2by3)';
doFH    = true; % true => FreemanHalton exact, false => chi-square
Tsummary2x3 = [Tsummary2x3; run2x3(myCat, catList, measure, doFH)];
doPlotExport(plot_option, save_dir, measure, 'svg');

%% 4) Animal ACC
myCat   = T1_epil_measures_upted.AnimalACCcat;
catList = {'Low','Mid','High'};
measure = 'Animal ACC (2by3)';
doFH    = true;
Tsummary2x3 = [Tsummary2x3; run2x3(myCat, catList, measure, doFH)];
doPlotExport(plot_option, save_dir, measure, 'svg');


%% 5) Animal RT
myCat   = T1_epil_measures_upted.AnimalRTcat;
catList = {'Fast','Moderate','Slow'};
measure = 'Animal RT (2by3)';
doFH    = true;  % e.g. chi-square
Tsummary2x3 = [Tsummary2x3; run2x3(myCat, catList, measure, doFH)];
doPlotExport(plot_option, save_dir, measure, 'svg');


%% 6) Symbol RT
myCat   = T1_epil_measures_upted.SymbolRTcat;
catList = {'Fast','Moderate','Slow'};
measure = 'Symbol RT (2by3)';
doFH    = true;
Tsummary2x3 = [Tsummary2x3; run2x3(myCat, catList, measure, doFH)];
doPlotExport(plot_option, save_dir, measure, 'svg');


%% 7) TLE side
myCat   = T1_epil_measures_upted.TLEside;
catList = {'Left','Right','Bilateral'};
measure = 'TLE side (2by3)';
doFH    = true;
Tsummary2x3 = [Tsummary2x3; run2x3(myCat, catList, measure, doFH)];
doPlotExport(plot_option, save_dir, measure, 'svg');


%% 8) EHQ
myCat   = T1_epil_measures_upted.EHQcat;
catList = {'Left','Right','Ambi'};
measure = 'EHQ (2by3)';
doFH    = true;
Tsummary2x3 = [Tsummary2x3; run2x3(myCat, catList, measure, doFH)];
doPlotExport(plot_option, save_dir, measure, 'svg');


%% 9) AED
myCat   = T1_epil_measures_upted.AEDcat;
catList = {'0','1','2plus'}; % example bins
measure = 'AED (2by3)';
doFH    = true;
Tsummary2x3 = [Tsummary2x3; run2x3(myCat, catList, measure, doFH)];
doPlotExport(plot_option, save_dir, measure, 'svg');


%% 10) LTGTC
myCat   = T1_epil_measures_upted.LTGTCcat;
catList = {'1-5','6-20','21plus'};
measure = 'LTGTC (2by3)';
doFH    = true;
Tsummary2x3 = [Tsummary2x3; run2x3(myCat, catList, measure, doFH)];
doPlotExport(plot_option, save_dir, measure, 'svg');


%% 11) SG
myCat   = T1_epil_measures_upted.SGcat;
catList = {'0','1-2','3plus'};
measure = 'SG freq (2by3)';
doFH    = true;
Tsummary2x3 = [Tsummary2x3; run2x3(myCat, catList, measure, doFH)];
doPlotExport(plot_option, save_dir, measure, 'svg');

%% 12) CP freq
myCat   = T1_epil_measures_upted.cp_freq_cat;
catList = {'6to10','1to5','0'};
measure = 'CP freq (2by3)';
doFH    = true;
Tsummary2x3 = [Tsummary2x3; run2x3(myCat, catList, measure, doFH)];
doPlotExport(plot_option, save_dir, measure, 'svg');

%% Finally, display
disp(Tsummary2x3);
writetable(Tsummary2x3,fullfile(save_dir, 'Fisher_Freeman_Halton_2by3.csv'));
