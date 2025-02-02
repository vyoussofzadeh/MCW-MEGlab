function [countsMatrix, pValsExact2x3, chi2P_vals, figHandle] = ecpfunc_plotDiscordantStackedBar_2x3(cfg)
% ECPFUNC_PLOTDISCORDANTSTACKEDBAR_2X3
%
% This function:
%  1) Takes a 'bestResultsTable' with columns: 'ROI' and 'discordColumn'.
%  2) Takes "myCategorical", an N×1 categorical with exactly 3 categories 
%     (e.g. 'Low','Mid','High').
%  3) For each ROI, we find "Discordant" vs "Concordant" => 2 groups
%     and tally how many are Low, Mid, High => a 2×3 table.
%  4) We do a stacked bar chart of # of Discordant subjects in {Low,Mid,High}.
%  5) If doFreemanHalton==true, we also run an exact 2×3 test (FreemanHalton).
%     If doFreemanHalton==false, we run a 2×3 Chi-square.
%
% USAGE:
%   cfg.bestResultsTable = bestResultsTable;
%   cfg.myCategorical    = T1_epil_measures_upted.AnimalACCcat; % or SymbolACCcat
%   cfg.discordColumn    = 'Gross_Discord_Subs';
%   cfg.categoryList     = {'Low','Mid','High'};
%   cfg.nSubjects        = 71 or 72, etc.
%   cfg.title            = 'Animal ACC 2x3';
%   cfg.doFreemanHalton  = true;  % if you want an exact p-value for each ROI
%
%   [counts, pVals2x3, chi2vals, figH] = ecpfunc_plotDiscordantStackedBar_2x3(cfg);
%
% OUTPUT:
%   countsMatrix    - [nROIs x 3] # of discordant subjects in each category
%   pValsExact2x3   - exact p-values from a custom 2×3 test (NaN if doFreemanHalton=false)
%   chi2P_vals      - chi-square p-values for 2×3 (NaN if doFreemanHalton=true)
%   figHandle       - handle to the stacked bar figure
%
% By [Your Name], [Date].

%% Unpack
bestResultsTable  = cfg.bestResultsTable;
myCategorical     = cfg.myCategorical;  % 3-level cat, e.g. 'Low','Mid','High'
discordColumn     = cfg.discordColumn;
categoryList      = cfg.categoryList;   % e.g. {'Low','Mid','High'}
nSubjects         = cfg.nSubjects;
if isfield(cfg,'title'), figTitle = cfg.title; else, figTitle="2x3 test"; end
if isfield(cfg,'doFreemanHalton'), doFH=cfg.doFreemanHalton; else, doFH=false; end

nCats  = length(categoryList);
allROI = bestResultsTable.ROI;
nROIs  = height(bestResultsTable);

% We'll store results
pValsExact2x3 = nan(nROIs,1);  % if doFreemanHalton==true
chi2P_vals    = nan(nROIs,1);  % if doFreemanHalton==false
countsMatrix  = zeros(nROIs, nCats);

allSubs = 1:nSubjects;

for i = 1:nROIs
    % (A) Discord vs Concord
    discordSubs  = bestResultsTable.(discordColumn){i};
    concordSubs  = setdiff(allSubs, discordSubs);

    % (B) Among the DISCORDANT group, how many are Low,Mid,High?
    cat_discord = myCategorical(discordSubs);
    for c = 1:nCats
        countsMatrix(i,c) = sum(cat_discord == categoryList{c});
    end

    % (C) Build 2×3 freq table: row1=Discord, row2=Concord; col=Low,Mid,High
    cat_concord = myCategorical(concordSubs);
    RowDiscord = zeros(1,nCats);
    RowConcord = zeros(1,nCats);
    for c = 1:nCats
        RowDiscord(c) = sum(cat_discord == categoryList{c});
        RowConcord(c) = sum(cat_concord == categoryList{c});
    end
    ContTable_2x3 = [RowDiscord; RowConcord];  % 2×3

    % (D) If doFreemanHalton => use exact 2x3 approach, else do Chi-square
    if doFH
        % Attempt a 2×3 exact test => user must provide a function e.g. freemanHalton2x3
        % or modify the code to use a general FreedmanHalton approach. We store the p-value:
        pValsExact2x3(i) = freemanHalton2x3(ContTable_2x3); 
        % ^ This is a placeholder function name. 
        %   You must implement or find user-submitted code for 2×3 exact test.
    else
        % Do a 2×3 Chi-square test
        chi2P_vals(i) = chiSquare2x3(ContTable_2x3);
    end
end

%% Plot stacked bar (Discordant group only)
figHandle = figure('Color','w','Name','Discordant 2×3','Position',[100 300 700 400]);
hBar = bar(countsMatrix, 'stacked');
colormapLines = lines(nCats);
for c = 1:min(nCats, size(colormapLines,1))
    hBar(c).FaceColor = colormapLines(c,:);
end
set(gca,'XTick',1:nROIs,'XTickLabel',allROI,'FontSize',10);
xlabel('ROI'); ylabel('# Discordant Subjects'); title(figTitle);
legend(categoryList,'Location','bestoutside');
box off;

%% Print results
for i = 1:nROIs
    fprintf('ROI: %s', allROI{i});
    if doFH && ~isnan(pValsExact2x3(i))
        fprintf(' | FreedmanHalton 2×3 p=%.4g', pValsExact2x3(i));
    elseif ~doFH && ~isnan(chi2P_vals(i))
        fprintf(' | Chi-square 2×3 p=%.4g', chi2P_vals(i));
    end
    fprintf('\n');
end

end

%% Subfunction: chiSquare2x3
function pVal = chiSquare2x3(ContTable_2x3)
% A simple 2×3 chi-square test of independence

observed = ContTable_2x3;  % 2×3
rowSum   = sum(observed,2);
colSum   = sum(observed,1);
N        = sum(observed,'all');
expected = (rowSum * colSum)/N;  % (2×1)*(1×3) => 2×3
chi2stat = sum( (observed - expected).^2 ./ expected,'all');
df       = (size(observed,1)-1)*(size(observed,2)-1); % => (2-1)*(3-1)=2
pVal     = 1-chi2cdf(chi2stat, df);
end
