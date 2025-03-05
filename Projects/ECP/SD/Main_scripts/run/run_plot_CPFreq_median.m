%% run_plot_CPFreq_median_fromCat.m
%
% 1) Convert your categorical "cp_freq_cat" to a numeric array "cpNumeric_72".
%    We'll map:
%        '0'       -> 0
%        '1to5'    -> 1
%        '6to10'   -> 2
%        '11plus'  -> 3
%        'Missing' -> NaN
% 2) Compute the median of cpNumeric_72 (ignoring NaNs).
% 3) Classify each subject's CP freq as <= median vs > median.
% 4) For each ROI, do a 2×2 Fisher test:
%      (Discordant vs Concordant) × (CP <= median vs CP > median).
% 5) Plot a 2-column stacked bar for the DISCORDANT group,
%    mentioning the median in the figure title.

%% Step 0: Convert cp_freq_cat -> numeric array cpNumeric_72

% Suppose T1_epil_measures.cp_freq_cat is a 72x1 categorical with categories:
%   '0','1to5','6to10','11plus' (maybe also 'Missing').
% We map them to 0,1,2,3,NaN for a median-based approach.

cpNumeric_72 = nan(72,1);  % pre-allocate

for i = 1:72
    switch char(cp_nums_72(i))  % 'cp_nums_72' is your categorical array
        case '0'
            cpNumeric_72(i) = 0;
        case '1to5'
            cpNumeric_72(i) = 1;
        case '6to10'
            cpNumeric_72(i) = 2;
        case '11plus'
            cpNumeric_72(i) = 3;
        case 'Missing'
            cpNumeric_72(i) = NaN;
        otherwise
            warning('Unrecognized CP freq category: %s', char(cp_nums_72(i)));
            cpNumeric_72(i) = NaN;
    end
end

%% A) Basic Setup
allROI   = bestResultsTable.ROI;     % cell array of ROI names
nROIs    = height(bestResultsTable); % number of rows
allSubs  = 1:72;                     % total # subjects [1..72]

% Which column for "discordant"?
discordColumn = 'Best_Discord_Subs'; 
% or 'Gross_Discord_Subs' if you have a "gross mismatch"

%% 1) Compute the median across cpNumeric_72
cp_median = nanmedian(cpNumeric_72);
fprintf('Median CP freq among 72 subjects = %.2f\n', cp_median);

%% B) Pre-allocate arrays for the stacked bar of DISCORDANT ONLY
% col1 = # with CP <= median, col2 = # with CP > median
n_leMed = zeros(nROIs,1);
n_gtMed = zeros(nROIs,1);

% We'll store p-values from the 2×2 Fisher test
fisherP_vals = nan(nROIs,1);

for i = 1:nROIs
    %% 1) Identify Discordant vs Concordant
    discordSubs = bestResultsTable.(discordColumn){i};
    concordSubs = setdiff(allSubs, discordSubs);

    %% 2) Subset cpNumeric_72 for these groups
    cp_discord = cpNumeric_72(discordSubs);
    cp_concord = cpNumeric_72(concordSubs);

    %% 3) For DISCORDANT: how many are <= median vs > median
    n_leMed(i) = sum(cp_discord <= cp_median);
    n_gtMed(i) = sum(cp_discord >  cp_median);

    %% 4) 2×2 Fisher: (Discord vs Concord) × (CP <= median vs CP > median)
    leMed_disc = sum(cp_discord <= cp_median);
    gtMed_disc = sum(cp_discord >  cp_median);
    leMed_conc = sum(cp_concord <= cp_median);
    gtMed_conc = sum(cp_concord >  cp_median);

    ContTable_2x2 = [leMed_disc, gtMed_disc;
                     leMed_conc, gtMed_conc];

    if sum(ContTable_2x2(1,:))>0 && sum(ContTable_2x2(2,:))>0
        [~, pVal] = fishertest(ContTable_2x2);
        fisherP_vals(i) = pVal;
    else
        fisherP_vals(i) = NaN;
    end

%     % Print
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('(Discord vs Concord) × (CP <= median vs CP > median):');
%     disp(ContTable_2x2);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
end

%% C) Plot the 2-column stacked bar for the DISCORDANT group
barData = [n_leMed, n_gtMed];  % Nx2

figure('Color','w','Name','Discordant By CP freq (<=median vs >median)', ...
       'Position',[100 300 700 400]);
hBar = bar(barData, 'stacked');

% Simple color scheme
colorSpec = [0.8 0.2 0.2;  0.2 0.6 0.8];  % red, teal
for c = 1:size(barData,2)
    hBar(c).FaceColor = colorSpec(c,:);
end

set(gca,'XTick',1:nROIs,'XTickLabel',allROI,'FontSize',9);
xlabel('ROI');
ylabel('# Discordant Patients');

% Mention the median in the figure title
medianTitle = sprintf('Discordant By CP Freq (<= median=%.2f vs > median)', cp_median);
title(medianTitle);

legend({'<= median','> median'}, 'Location','bestoutside');
box off;

%% D) Print final p-values
for i = 1:nROIs
    fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
end
