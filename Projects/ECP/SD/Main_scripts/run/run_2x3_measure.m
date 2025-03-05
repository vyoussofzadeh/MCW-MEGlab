function Tloc = run_2x3_measure(bestResultsTable, myCat, catList, measureName, doFH)
% RUN_2X3_MEASURE: helper to reduce repeated code
%   - builds the cfg
%   - calls ecpfunc_plotDiscordantStackedBar_2x3
%   - returns Nx4 table => (Measure, Comparison, ROI, pVal)

% cfg = [];
% cfg.bestResultsTable = bestResultsTable;
% cfg.myCategorical    = myCat;
% cfg.discordColumn    = 'Gross_Discord_Subs';  % or 'Best_Discord_Subs'
% cfg.categoryList     = catList;
% cfg.nSubjects        = length(myCat); % assuming same # of subjects
% cfg.title            = measureName;
% cfg.doFreemanHalton  = doFH;
% 
% [counts2x3, pVals2x3, chi2vals, ~] = ecpfunc_plotDiscordantStackedBar_2x3(cfg);

cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical    = myCat;
cfg.discordColumn    = 'Gross_Discord_Subs';  % or whichever column has the discord list
cfg.categoryList     = catList;
cfg.nSubjects        = length(myCat);  % or length(T1_epil_measures_upted.SubjectID)
cfg.title            = measureName;
cfg.doFreemanHalton  = doFH;  % if true => FreemanHalton, else => 2×3 Chi-square

[counts2x3, pVals2x3, chi2vals, ~] = ecpfunc_plotDiscordantStackedBar_2x3(cfg);

nROIs = height(bestResultsTable);
% Tloc = table( ...
%     repmat(measureName,nROIs,1), ...
%     repmat(catList,nROIs,1), ...
%     bestResultsTable.ROI, ...
%     pVals2x3, ...
%     'VariableNames',{'Measure','Comparison','ROI','pVal'});

comparisonStr = sprintf('2by3(Discord vs Concord), classes={%s}', ...
    strjoin(catList, ','));

Tloc = table( ...
    repmat(measureName,nROIs,1), ...
    repmat(comparisonStr,nROIs,1), ...
    bestResultsTable.ROI, ...
    pVals2x3, ...
    'VariableNames',{'Measure','Comparison','ROI','pVal'});
end
