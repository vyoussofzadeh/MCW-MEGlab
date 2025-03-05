%% run_plot_SGfreq_median_fromCat.m
%
% 1) Convert your categorical sg_nums_72 to a numeric array sgNumeric_72
%    mapping e.g. '0'->0, '1to2'->1, '3plus'->2, etc.
% 2) Compute the median (ignoring NaN).
% 3) Classify each subject as <= median vs > median.
% 4) For each ROI, do a 2×2 Fisher test: (Discord vs Concord) × (SG <= median vs SG > median).
% 5) Plot a 2-column stacked bar for the DISCORDANT group, and mention the median in the figure title.

%% Step 0: Convert sg_nums_72 (categorical) -> numeric array
% Suppose sg_nums_72 is a 72x1 categorical with categories like:
% {'0','1to2','3plus'} and possibly 'Missing', etc.
% We'll map them to numeric:
%   '0'     -> 0
%   '1to2'  -> 1
%   '3plus' -> 2
%   'Missing' -> NaN
% You can adjust or expand if you have more categories.

sgNumeric_72 = nan(72,1);  % pre-allocate

for i = 1:72
    switch char(sg_nums_72(i))
        case '0'
            sgNumeric_72(i) = 0;
        case '1to2'
            sgNumeric_72(i) = 1;
        case '3plus'
            sgNumeric_72(i) = 2;
        case 'Missing'
            sgNumeric_72(i) = NaN;
        otherwise
            warning('Unrecognized SG category: %s', char(sg_nums_72(i)));
            sgNumeric_72(i) = NaN;
    end
end

%% A) Basic checks
allROI = bestResultsTable.ROI;      % cell array of ROI names
nROIs  = height(bestResultsTable);  % number of rows
allSubs = 1:72;                     % total # subjects (1..72)

% Which column for "discordant"?
discordColumn = 'Best_Discord_Subs';
% or 'Gross_Discord_Subs' if you want a "gross mismatch"

%% 1) Compute the median of the numeric array (ignoring NaNs)
sg_median = nanmedian(sgNumeric_72);
fprintf('Median SG frequency among 72 subjects = %.2f\n', sg_median);

%% 2) We'll store 2 columns for the stacked bar in the DISCORDANT group:
% col1 = # with SG <= median, col2 = # with SG > median
n_leMed = zeros(nROIs,1);
n_gtMed = zeros(nROIs,1);

% We'll store p-values from the 2×2 Fisher test
fisherP_vals = nan(nROIs,1);

for i = 1:nROIs
    %% A) Identify Discordant vs. Concordant
    discordSubs = bestResultsTable.(discordColumn){i};
    concordSubs = setdiff(allSubs, discordSubs);

    %% B) Subset sgNumeric_72 for these groups
    sg_discord = sgNumeric_72(discordSubs);
    sg_concord = sgNumeric_72(concordSubs);

    %% C) For DISCORDANT: how many are <= median vs > median
    n_leMed(i) = sum(sg_discord <= sg_median);
    n_gtMed(i) = sum(sg_discord >  sg_median);

    %% 2×2 table: (Discord vs Concord) x (SG <= median vs SG > median)
    leMed_disc = sum(sg_discord <= sg_median);
    gtMed_disc = sum(sg_discord >  sg_median);

    leMed_conc = sum(sg_concord <= sg_median);
    gtMed_conc = sum(sg_concord >  sg_median);

    ContTable_2x2 = [leMed_disc, gtMed_disc;
                     leMed_conc, gtMed_conc];

    % Only run fisher if each row > 0 total
    if sum(ContTable_2x2(1,:))>0 && sum(ContTable_2x2(2,:))>0
        [~, pVal] = fishertest(ContTable_2x2);
        fisherP_vals(i) = pVal;
    else
        fisherP_vals(i) = NaN;
    end

%     % Print table & p-value
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('(Discord vs Concord) × (SG <= median vs SG > median):');
%     disp(ContTable_2x2);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
end

%% 3) Plot the 2-column stacked bar for the DISCORDANT group
barData = [n_leMed, n_gtMed];  % Nx2

figure('Color','w','Name','Discordant By SG freq (<=median vs >median)', ...
    'Position',[100 300 700 400]);
hBar = bar(barData, 'stacked');

% Simple color scheme
colorSpec = [0.8 0.2 0.2; 0.2 0.6 0.8];  % red, teal
for c = 1:size(barData,2)
    hBar(c).FaceColor = colorSpec(c,:);
end

set(gca,'XTick',1:nROIs,'XTickLabel',allROI,'FontSize',9);
xlabel('ROI');
ylabel('# Discordant Patients');

% Mention the median in the figure title:
medianTitle = sprintf('Discordant by SG Frequency (<= median=%.2f vs > median)', sg_median);
title(medianTitle);

legend({'<= median','> median'}, 'Location','bestoutside');
box off;

%% 4) Print final p-values
for i = 1:nROIs
    fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
end
