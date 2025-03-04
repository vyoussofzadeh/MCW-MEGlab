function [countsMatrix, fisherP_vals, oddsRatio_vals, chi2P_vals, figHandle] = ecpfunc_plotDiscordantStackedBar(cfg)
% ECPFUNC_PLOTDISCORDANTSTACKEDBAR
%
% This function:
%  1) Takes a table "bestResultsTable" with columns: 'ROI' and 'discordColumn'.
%  2) Takes "myCategorical", an N×1 categorical array for your N subjects
%     (like T1_epil_measures_upted.SymbolACCcat). The default was 72, but
%     now you can have "nSubjects" = 71 or any integer.
%  3) "categoryList" is a cell array of categories to display in the stacked bar.
%     If empty ([]), we use categories(myCategorical).
%  4) If "doFisher" == true, we do a custom 2×2 test by collapsing categories
%     into "binA" vs. "binB". e.g. binA = {'Low','Mid'}, binB = {'High'}.
%  5) We create a stacked bar chart of how many DISCORDANT subjects in each category,
%     one bar per ROI, and (optionally) return p-values & odds ratios for the 2×2 test.
%  6) "nSubjects" is the total # of subjects your code references (could be 71, 72, etc.).
%
% USAGE:
%   [counts, pVals, ORvals, chi2P_vals, hFig] = ecpfunc_plotDiscordantStackedBar( ...
%       bestResultsTable, ...
%       T1_epil_measures_upted.SymbolACCcat, ...
%       'Gross_Discord_Subs', ...
%       {'Low','Mid','High'}, ...
%       true, ...
%       {'Low','Mid'}, ...      % binA
%       {'High'}, ...           % binB
%       71);
%
% Where we do a 2×2 test by lumping 'Low','Mid' into binA, 'High' into binB,
% ignoring categories not in binA or binB for the Fisher and Chi-Square tests.
%
% INPUTS:
%   bestResultsTable - table with 'ROI' (nROIs x 1) and e.g. 'Gross_Discord_Subs' columns
%   myCategorical    - nSubjects×1 categorical array for your dataset (N could be 71 or 72, etc.)
%   discordColumn    - string, the name of the column in bestResultsTable that contains the discord subs
%   categoryList     - cell array of categories to display in the stacked bar
%                      (if empty, we use categories(myCategorical)).
%   doFisher         - boolean; if true, do a 2×2 fisher test
%   binA             - cell array of category names that define one bin for the fisher test
%   binB             - cell array of category names that define the other bin for the fisher test
%   nSubjects        - integer specifying how many total subjects you have (could be 71, 72, etc.)
%
% OUTPUTS:
%   countsMatrix     - [nROIs x nCats], # of discordant subjects in each category
%   fisherP_vals     - [nROIs x 1], p-values from the 2×2 fisher test (NaN if doFisher=false)
%   oddsRatio_vals   - [nROIs x 1], odds ratio from fisher (NaN if doFisher=false)
%   chi2P_vals       - [nROIs x 1], p-values from the 2×2 chi-square test (NaN if doFisher=false)
%   figHandle        - handle to the created figure
%
% By [Your Name], [Date].

bestResultsTable = cfg.bestResultsTable;
myCategorical = cfg.myCategorical;
discordColumn = cfg.discordColumn;
categoryList = cfg.categoryList;
doFisher = cfg.doFisher;
binA = cfg.binA;
binB  =cfg.binB;
nSubjects = cfg.nSubjects;

nCats = length(categoryList);
allROI = bestResultsTable.ROI;
nROIs  = height(bestResultsTable);

% We'll store p-values and OR
fisherP_vals   = nan(nROIs,1);
oddsRatio_vals = nan(nROIs,1);
chi2P_vals     = nan(nROIs,1); % New output for chi-square p-values

% We'll store counts for the stacked bar
countsMatrix = zeros(nROIs, nCats);

allSubs = 1:nSubjects;  % total subjects

for i = 1:nROIs
    %% 1) Get discordant subset for this ROI
    discordSubs = bestResultsTable.(discordColumn){i};
    % Concordant
    concordSubs = setdiff(allSubs, discordSubs);

    %% 2) For the DISCORDANT group, retrieve the categories
    cat_discord = myCategorical(discordSubs);

    %% 3) Count how many fall into each category
    for c = 1:nCats
        countsMatrix(i,c) = sum(cat_discord == categoryList{c});
    end

    %% 4) If doFisher is true, do a custom 2×2 test
    if doFisher && ~isempty(binA) && ~isempty(binB)
        cat_concord = myCategorical(concordSubs);

        % We'll sum up binA and binB among discord & concord
        nA_disc = sum(ismember(cat_discord, binA));
        nB_disc = sum(ismember(cat_discord, binB));

        nA_conc = sum(ismember(cat_concord, binA));
        nB_conc = sum(ismember(cat_concord, binB));

        ContTable_2x2 = [nA_disc, nB_disc;
                         nA_conc, nB_conc];

        % If each row has >0 total
        if all(sum(ContTable_2x2,2) > 0)
            % Fisher's exact test
            [~, pVal, stats] = fishertest(ContTable_2x2);
            fisherP_vals(i)   = pVal;
            oddsRatio_vals(i) = stats.OddsRatio;

            % Chi-square test
            [~, chi2P_vals(i)] = chi2test(ContTable_2x2); % Call chi2test function
        else
            fisherP_vals(i)   = NaN;
            oddsRatio_vals(i) = NaN;
            chi2P_vals(i)     = NaN;
        end
    end
end

%% Now create the stacked bar of # DISCORDANT subjects by category
figHandle = figure('Color','w','Name','Discordant Stacked Bar','Position',[100 300 700 400]);
hBar = bar(countsMatrix, 'stacked');

% Assign some colors if you'd like
colormapLines = lines(nCats); % or define your own color map
for c = 1:min(nCats, size(colormapLines,1))
    hBar(c).FaceColor = colormapLines(c,:);
end

% Label x-axis with ROI
set(gca,'XTick',1:nROIs,'XTickLabel',allROI,'FontSize',10);
xlabel('ROI');
ylabel('# Discordant Subjects');
title(cfg.title)
% title('Discordant Subjects by Category');
legend(categoryList, 'Location','bestoutside');
box off;

%% (Optional) Print fisher and chi-square results
for i = 1:nROIs
%     fprintf('ROI: %s', allROI{i});
    if doFisher && ~isnan(fisherP_vals(i))
%         fprintf(' | Fisher p=%.4g, OR=%.3f', fisherP_vals(i), oddsRatio_vals(i));
    end
    if doFisher && ~isnan(chi2P_vals(i))
%         fprintf(' | Chi-Square p=%.4g', chi2P_vals(i));
    end
%     fprintf('\n');
end

end

%% Helper function for chi-square test
function [chi2stat, pVal] = chi2test(observed)
    % chi2test - Performs a chi-square test of independence.
    %
    % Input:
    %   observed - A 2×2 contingency table of observed frequencies.
    %
    % Output:
    %   chi2stat - The chi-square test statistic.
    %   pVal     - The p-value of the test.

    % Calculate the expected frequencies
    rowTotals = sum(observed, 2);
    colTotals = sum(observed, 1);
    total = sum(rowTotals);
    expected = (rowTotals * colTotals) / total;

    % Calculate the chi-square statistic
    chi2stat = sum((observed - expected).^2 ./ expected, 'all');

    % Calculate degrees of freedom
    df = 1; % For a 2×2 table, df = (2-1)*(2-1) = 1

    % Calculate the p-value
    pVal = 1 - chi2cdf(chi2stat, df);
end