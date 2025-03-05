%% run_plot_LTGTC_median_fromCat.m
%
% 1) Convert your categorical ltgtc_nums_72 -> numeric array ltgtcNumeric_72
%    mapping '0','1-5','6-20','21plus' -> 0,1,2,3 and 'Missing' -> NaN.
% 2) Compute the median of ltgtcNumeric_72 (ignoring NaNs).
% 3) Classify each subject's numeric LTGTC as <= median vs > median.
% 4) For each ROI, do a 2×2 Fisher test:
%      (Discordant vs Concordant) × (LTGTC <= median vs LTGTC > median).
% 5) Plot a 2-column stacked bar for the DISCORDANT group, mentioning the median value.

%% Step 0: Convert ltgtc_nums_72 (categorical) -> numeric

ltgtcNumeric_72 = nan(72,1);  % pre-allocate

for i = 1:72
    switch char(ltgtc_nums_72(i))
        case '0'
            ltgtcNumeric_72(i) = 0;
        case '1-5'
            ltgtcNumeric_72(i) = 1;
        case '6-20'
            ltgtcNumeric_72(i) = 2;
        case '21plus'
            ltgtcNumeric_72(i) = 3;
        case 'Missing'
            ltgtcNumeric_72(i) = NaN;
        otherwise
            warning('Unrecognized LTGTC category: %s', char(ltgtc_nums_72(i)));
            ltgtcNumeric_72(i) = NaN;
    end
end

%% (1) Basic Setup
allROI   = bestResultsTable.ROI;   % cell array of ROI names
nROIs    = height(bestResultsTable);
allSubs  = 1:72;                  % total # of subjects in your main set

% Which column for "discordant"?
discordColumn = 'Best_Discord_Subs';  
% or 'Gross_Discord_Subs' if you want a "gross mismatch"

%% (2) Compute the median across the 72 subjects (ignoring NaNs)
ltgtc_median = nanmedian(ltgtcNumeric_72);
fprintf('Median numeric LTGTC among 72 subjects = %.1f\n', ltgtc_median);

%% (3) We'll store 2 columns for the stacked bar (Discordant only):
% col1 => # with LTGTC <= median
% col2 => # with LTGTC >  median
n_leMed = zeros(nROIs,1);
n_gtMed = zeros(nROIs,1);

% We'll store p-values from the 2×2 Fisher test
fisherP_vals = nan(nROIs,1);

for i = 1:nROIs
    % (A) Identify Discordant vs. Concordant for this ROI
    discordSubs = bestResultsTable.(discordColumn){i};
    concordSubs = setdiff(allSubs, discordSubs);

    % (B) Subset ltgtcNumeric_72 for these groups
    ltgtc_discord = ltgtcNumeric_72(discordSubs);
    ltgtc_concord = ltgtcNumeric_72(concordSubs);

    % (C) Count how many are <= median vs. > median in DISCORDANT group
    n_leMed(i) = sum(ltgtc_discord <= ltgtc_median);
    n_gtMed(i) = sum(ltgtc_discord >  ltgtc_median);

    % Build the 2×2 table:
    leMed_disc = sum(ltgtc_discord <= ltgtc_median);
    gtMed_disc = sum(ltgtc_discord >  ltgtc_median);

    leMed_conc = sum(ltgtc_concord <= ltgtc_median);
    gtMed_conc = sum(ltgtc_concord >  ltgtc_median);

    ContTable_2x2 = [leMed_disc, gtMed_disc;
                     leMed_conc, gtMed_conc];

    if sum(ContTable_2x2(1,:))>0 && sum(ContTable_2x2(2,:))>0
        [~, pVal] = fishertest(ContTable_2x2);
        fisherP_vals(i) = pVal;
    else
        fisherP_vals(i) = NaN;
    end

%     % Print table & p-value
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('(Discord vs Concord) × (LTGTC <= median vs LTGTC > median):');
%     disp(ContTable_2x2);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
end

%% (4) Plot the 2-column stacked bar for the DISCORDANT group
barData = [n_leMed, n_gtMed];  % Nx2

figure('Color','w','Name','Discordant By LTGTC (<=median vs >median)', ...
       'Position',[100 300 700 400]);
hBar = bar(barData, 'stacked');

% Color scheme
colorSpec = [0.2 0.6 0.8; 0.8 0.2 0.2];  % teal, red
for c = 1:2
    hBar(c).FaceColor = colorSpec(c,:);
end

set(gca,'XTick',1:nROIs,'XTickLabel',allROI,'FontSize',10);
xlabel('ROI');
ylabel('# Discordant Patients');

% Mention the median in the figure title
medianTitle = sprintf('Discordant by LTGTC (<= median=%.1f vs > median)', ltgtc_median);
title(medianTitle);

legend({'<= median','> median'}, 'Location','bestoutside');
box off;

%% (5) Print final p-values
for i = 1:nROIs
    fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
end

