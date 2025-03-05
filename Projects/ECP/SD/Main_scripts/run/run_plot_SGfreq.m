%% run_plot_SGfreq_simplified.m
%
% This script merges the "SGcat" categories into two bins:
%   "0"   vs.   "=1"  (where "=1" includes '1-2','3plus', etc.)
%
% Then we run a 2×2 Fisher test: (Discordant vs Concordant) × (SG=0 vs SG=1).
% Finally, we plot a 2-column stacked bar for the DISCORDANT group.

%% A) Basic checks
allROI = bestResultsTable.ROI;        % cell array of ROI names
nROIs  = height(bestResultsTable);    % # rows in bestResultsTable
allSubs = 1:72;                       % total # subjects (adjust as needed)

% Decide which discord column to use
discordColumn = 'Best_Discord_Subs';  
% or 'Gross_Discord_Subs' if you want "gross mismatch"

%% B) Pre-allocate arrays for the stacked bar of DISCORDANT ONLY
% col1 = # with SG=0, col2 = # with SG=1
n_zero = zeros(nROIs,1);
n_plus = zeros(nROIs,1);

% We'll store p-values from the 2×2 Fisher test
fisherP_vals = nan(nROIs,1);

for i = 1:nROIs
    % 1) Identify discordant vs. concordant
    discordSubs = bestResultsTable.(discordColumn){i};
    concordSubs = setdiff(allSubs, discordSubs);

    % 2) Subset T1_epil_measures_upted.SGcat for these groups
    sg_discord = T1_epil_measures_upted.SGcat(discordSubs);
    sg_concord = T1_epil_measures_upted.SGcat(concordSubs);

    % 3) For DISCORDANT: how many are "0" vs. "=1" (i.e. anything else)
    n_zero(i) = sum(sg_discord == "0");
    n_plus(i) = sum(sg_discord ~= "0");  % lumps '1-2','3plus', etc.

    % ----- 2×2 Fisher: (Discord vs Concord) × (SG=0 vs SG=1) -----
    zero_disc   = sum(sg_discord == "0");
    plus_disc   = sum(sg_discord ~= "0");
    zero_conc   = sum(sg_concord == "0");
    plus_conc   = sum(sg_concord ~= "0");

    ContTable_2x2 = [zero_disc, plus_disc;
                     zero_conc, plus_conc];

    % Only run fisher if both rows have > 0 total
    if sum(ContTable_2x2(1,:))>0 && sum(ContTable_2x2(2,:))>0
        [~, pVal] = fishertest(ContTable_2x2);
        fisherP_vals(i) = pVal;
    else
        fisherP_vals(i) = NaN;
    end

    % Print
    fprintf('\nROI #%d = %s\n', i, allROI{i});
    disp('(Discord vs Concord) × (SG=0 vs SG=1) :');
    disp(ContTable_2x2);
    fprintf('p-value = %.4g\n', fisherP_vals(i));
end

%% C) Plot the 2-column stacked bar for DISCORDANT ONLY
barData = [n_zero, n_plus];   % Nx2

figure('Color','w','Name','Discordant By SGfreq (0 vs >=1)','Position',[100 300 700 400]);
hBar = bar(barData, 'stacked');

% Simple color scheme
colorSpec = [0.8 0.2 0.2;  0.2 0.6 0.8]; % e.g. red, teal
for c = 1:size(barData,2)
    hBar(c).FaceColor = colorSpec(c,:);
end

set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize',9);
xlabel('ROI');
ylabel('# Discordant Patients');
title('Discordant Subjects by SG Frequency (0 vs. >=1)');
legend({'0','=1'}, 'Location','bestoutside');
box off;

%% Print final p-values
for i = 1:nROIs
    fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
end


% %% run_plot_SGfreq.m
% % Plots how many DISCORDANT subjects fall into each category of "SGcat"
% % (e.g. '0','1-2','3plus') for each ROI, and does a 2×n Fisher test
% % of (Discordant vs Concordant) × {SGcat categories}.
% 
% % A) Basic checks
% allROI = bestResultsTable.ROI;     % cell array of ROI names
% nROIs  = height(bestResultsTable); % number of ROIs
% 
% % Suppose you have 72 total subjects
% allSubs = 1:72;  
% 
% % B) Grab the categories from T1_epil_measures_upted.SGcat
% SGcats = categories(T1_epil_measures_upted.SGcat);  % e.g. {'0','1to2','3plus'}
% nCats  = length(SGcats);
% 
% % C) Pre-allocate a [nROIs x nCats] matrix for the DISCORDANT group
% SG_counts = zeros(nROIs, nCats);
% 
% % (Optional) store p-values for the multi-way Fisher test
% fisherP_vals = nan(nROIs,1);
% 
% for i = 1:nROIs
%     % 1) Identify discordant vs. concordant for this ROI
% %     discordSubs = bestResultsTable.Best_Discord_Subs{i};
%     discordSubs = bestResultsTable.Gross_Discord_Subs{i};
% 
%     concordSubs = setdiff(allSubs, discordSubs);
% 
%     % 2) Extract SGcat for these DISCORDANT subjects
%     sgVals_discord = T1_epil_measures_upted.SGcat(discordSubs);
% 
%     % 3) Count how many fall into each SGcat category among discordant
%     for c = 1:nCats
%         SG_counts(i,c) = sum(sgVals_discord == SGcats{c});
%     end
% 
%     % 4) Build a 2×n table for Fisher's exact if each row > 0
%     sgVals_concord = T1_epil_measures_upted.SGcat(concordSubs);
% 
%     ContTable = zeros(2, nCats);
%     for c = 1:nCats
%         ContTable(1,c) = sum(sgVals_discord  == SGcats{c});
%         ContTable(2,c) = sum(sgVals_concord  == SGcats{c});
%     end
% 
%     if sum(ContTable(1,:))>0 && sum(ContTable(2,:))>0
%         try
%             [~, pVal] = fishertest(ContTable, 'Tail','both');
%             fisherP_vals(i) = pVal;
%         catch ME
%             warning('Fisher test error in ROI %s: %s', allROI{i}, ME.message);
%             fisherP_vals(i) = NaN;
%         end
%     end
% 
%     % Print results
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('2×n Fisher test (Discordant/Concordant × SG categories):');
%     disp(ContTable);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
% end
% 
% %% D) Plot stacked bar for DISCORDANT
% figure('Color','w','Name','Discordant By SGfreq','Position',[100 300 700 400]);
% hBar = bar(SG_counts, 'stacked');
% 
% % E) Optional color scheme
% colorMap = [ ...
%     0.8 0.2 0.2;  ... red (if you have 1st category '0')
%     0.5 0.5 0.5;  ... gray (e.g. '1-2')
%     0.2 0.6 0.8;  ... teal (e.g. '3plus')
%     0.4 0.2 0.6;  ... if you have a 4th category
% ];
% 
% for c = 1:min(nCats,size(colorMap,1))
%     hBar(c).FaceColor = colorMap(c,:);
% end
% 
% % F) Axes, labels, legend
% set(gca, 'XTick',1:nROIs, 'XTickLabel', allROI, 'FontSize',9);
% xlabel('ROI');
% ylabel('# Discordant Patients');
% title('Discordant Subjects by SG Frequency (SGcat)');
% legend(SGcats, 'Location','bestoutside');
% box off;
% 
% %% Print final p-values
% for i = 1:nROIs
%     fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
% end
