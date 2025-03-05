%% run_plot_AED_median_fromCat.m
%
% 1) Convert your categorical aedCount_72 -> numeric array aedCountNumeric_72
%    mapping '0','1','2','3plus' -> 0,1,2,3.
% 2) Compute the median of aedCountNumeric_72.
% 3) Classify each subject's AED count as <= median vs > median.
% 4) For each ROI, do a 2×2 Fisher test:
%      (Discordant vs Concordant) × (AED <= median vs AED > median)
% 5) Plot a 2-column stacked bar for the DISCORDANT group
%    and include the median value in the figure title.

%% Step 0: Convert aedCount_72 (categorical) to numeric
aedCountNumeric_72 = nan(72,1);  % pre-allocate

for i = 1:72
    switch char(aedCount_72(i))
        case '0'
            aedCountNumeric_72(i) = 0;
        case '1'
            aedCountNumeric_72(i) = 1;
        case '2'
            aedCountNumeric_72(i) = 2;
        case '3plus'
            aedCountNumeric_72(i) = 3;
        otherwise
            warning('Unrecognized category: %s', char(aedCount_72(i)));
            aedCountNumeric_72(i) = NaN;
    end
end

%% 1) Basic Setup
allROI   = bestResultsTable.ROI;     
nROIs    = height(bestResultsTable); 
allSubs  = 1:72;                    
discordColumn = 'Best_Discord_Subs'; 
% or 'Gross_Discord_Subs' if you prefer "gross mismatch"

%% 2) Compute the median across the 72 subjects
medAED = nanmedian(aedCountNumeric_72);
fprintf('Median (numeric) AEDCount among the 72 patients = %.1f\n', medAED);

%% 3) We'll store 2 columns for the stacked bar (Discordant only):
n_leMed = zeros(nROIs,1);  % # of subjects with AED <= median
n_gtMed = zeros(nROIs,1);  % # of subjects with AED > median

% We'll store p-values from the 2×2 Fisher test
fisherP_vals = nan(nROIs,1);

for i = 1:nROIs
    %% A) Identify Discordant vs. Concordant for this ROI
    discordSubs = bestResultsTable.(discordColumn){i};
    concordSubs = setdiff(allSubs, discordSubs);

    %% B) Extract numeric AED counts for these subsets
    aed_discord = aedCountNumeric_72(discordSubs);
    aed_concord = aedCountNumeric_72(concordSubs);

    %% C) For the DISCORDANT group, count <= median vs. > median
    n_leMed(i) = sum(aed_discord <= medAED);
    n_gtMed(i) = sum(aed_discord >  medAED);

    %% 2×2 table for Fisher:
    leMed_disc = sum(aed_discord <= medAED);
    gtMed_disc = sum(aed_discord >  medAED);

    leMed_conc = sum(aed_concord <= medAED);
    gtMed_conc = sum(aed_concord >  medAED);

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
%     disp('(Discord vs Concord) × (AED <= median vs AED > median):');
%     disp(ContTable_2x2);
%     fprintf('p-value = %.4g\n', pVal);
end

%% 4) Plot the 2-column stacked bar for the DISCORDANT group
barData = [n_leMed, n_gtMed];  % Nx2

figure('Color','w','Name','Discordant By AED (<=median vs >median)',...
       'Position',[100 300 700 400]);
hBar = bar(barData, 'stacked');

% Color scheme
colorSpec = [0.8 0.2 0.2;  0.2 0.6 0.8];  % red, teal
for c = 1:2
    hBar(c).FaceColor = colorSpec(c,:);
end

set(gca,'XTick',1:nROIs,'XTickLabel',allROI,'FontSize',10);
xlabel('ROI');
ylabel('# Discordant Patients');

% Mention the median value in the figure title:
medianTitle = sprintf('Discordant Subjects by AED Count (<= median=%.1f vs > median)', medAED);
title(medianTitle);

legend({'<= median','> median'}, 'Location','bestoutside');
box off;

%% 5) Print final p-values
for i = 1:nROIs
    fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
end

