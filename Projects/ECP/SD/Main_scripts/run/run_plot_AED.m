%% run_plot_AED_0to1_vs_over1.m
%
% Merges the AED categories into two bins:
%   "0-1"  vs.  ">1"
% Then we do a 2×2 Fisher test:
%   (Discordant vs Concordant) × (AED=01 vs. AED>1)
% Finally, a 2-column stacked bar for the Discordant group only.

%% 1) Basic Setup
allROI = bestResultsTable.ROI;      % cell array of ROI names
nROIs  = height(bestResultsTable);  % number of rows in bestResultsTable
allSubs = 1:72;                     % total # subjects (adjust if needed)

% Which column do you use for "discordant"?
discordColumn = 'Best_Discord_Subs';
% or 'Gross_Discord_Subs' if you prefer "gross mismatch"

%% 2) We'll store 2 columns for the stacked bar in the Discordant group:
% col1 = "01", col2 = ">1"
n_0to1 = zeros(nROIs,1);
n_over1= zeros(nROIs,1);

% We'll store p-values from the 2×2 Fisher test
fisherP_vals = nan(nROIs,1);

for i = 1:nROIs
    %% A) Identify Discordant vs. Concordant
    discordSubs = bestResultsTable.(discordColumn){i};
    concordSubs = setdiff(allSubs, discordSubs);

    %% B) Subset T1_epil_measures.AEDcat for these groups
    aed_discord = T1_epil_measures_upted.AEDcat(discordSubs);
    aed_concord = T1_epil_measures_upted.AEDcat(concordSubs);

    %% C) For the DISCORDANT group, count how many are "01" vs. ">1"
    % "01" means aed_discord is '0' or '1'
    % ">1" means aed_discord is '2' or '3plus'
    n_0to1(i) = sum(aed_discord=="0" | aed_discord=="1");
    n_over1(i)= sum(aed_discord=="2" | aed_discord=="3plus");

    %% D) 2×2 Fisher: (Discord vs Concord) × (AED=01 vs AED>1)
    zeroOne_disc = sum(aed_discord=="0" | aed_discord=="1");
    over1_disc   = sum(aed_discord=="2" | aed_discord=="3plus");

    zeroOne_conc = sum(aed_concord=="0" | aed_concord=="1");
    over1_conc   = sum(aed_concord=="2" | aed_concord=="3plus");

    ContTable_2x2 = [zeroOne_disc,  over1_disc;
                     zeroOne_conc, over1_conc];

    % Only run fisher if each row > 0 total
    if sum(ContTable_2x2(1,:))>0 && sum(ContTable_2x2(2,:))>0
        [~, pVal] = fishertest(ContTable_2x2);
        fisherP_vals(i) = pVal;
    else
        fisherP_vals(i) = NaN;
    end

%     % Print
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('(Discord vs Concord) × (AED=0-1 vs AED>1):');
%     disp(ContTable_2x2);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
end

%% 3) Plot the 2-column stacked bar (Discordant only)
barData = [n_0to1, n_over1];  % Nx2

figure('Color','w','Name','Discordant By AED (0-1 vs >1)','Position',[100 300 700 400]);
hBar = bar(barData, 'stacked');

% Optional color scheme for 2 columns
colorSpec = [0.8 0.2 0.2;  0.2 0.6 0.8];  % red, teal
for c = 1:size(barData,2)
    hBar(c).FaceColor = colorSpec(c,:);
end

set(gca,'XTick',1:nROIs,'XTickLabel',allROI,'FontSize',10);
xlabel('ROI');
ylabel('# Discordant Patients');
title('Discordant Subjects by AED Count (0-1 vs >1)');
legend({'0-1','>1'}, 'Location','bestoutside');
box off;

%% 4) Print final p-values
for i = 1:nROIs
    fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
end

% allROI = bestResultsTable.ROI;       % cell array of ROI names
% nROIs  = height(bestResultsTable);   % number of rows (ROIs)
% 
% allSubs = 1:72;  % total # subjects (adjust if needed)
% 
% % We store the stacked bar counts for DISCORDANT only
% n0      = zeros(nROIs,1);  % # of subjects with '0' AED
% n1      = zeros(nROIs,1);  % # of subjects with '1'
% n2      = zeros(nROIs,1);  % # of subjects with '2'
% n3plus  = zeros(nROIs,1);  % # of subjects with '3plus'
% 
% % We'll store p-values from 2x2 Fisher test (discord vs. conc) x (0 vs. >=1)
% fisherP_vals = nan(nROIs,1);
% 
% for i = 1:nROIs
%     % 1) Grab the "discordant" subset for this ROI.
%     %    Or "grossly" discordant if you changed that above:
%     %    discordSubs = bestResultsTable.Gross_Discord_Subs{i};
%     discordSubs  = bestResultsTable.Best_Discord_Subs{i};
% 
%     % 2) The remainder are "concordant"
%     concordSubs  = setdiff(allSubs, discordSubs);
% 
%     % 3) Subset T1_epil_measures_upted.AEDcat for these groups
%     aed_discord  = T1_epil_measures_upted.AEDcat(discordSubs);
%     aed_concord  = T1_epil_measures_upted.AEDcat(concordSubs);
% 
%     % 4) For the stacked bar, we count only the DISCORDANT group
%     n0(i)     = sum(aed_discord == "0");
%     n1(i)     = sum(aed_discord == "1");
%     n2(i)     = sum(aed_discord == "2");
%     n3plus(i) = sum(aed_discord == "3plus");
% 
%     % ------- 2×2 Fisher: "0" vs. "=1" across Discordant vs Concordant ------
%     % Row 1 = Discordant, row 2 = Concordant
%     % Col 1 = # subjects with '0', col 2 = # subjects with '=1'
% 
%     n0_disc   = sum(aed_discord == "0");
%     nNon0_disc= sum(aed_discord ~= "0");   % i.e. "1","2","3plus"
%     n0_conc   = sum(aed_concord == "0");
%     nNon0_conc= sum(aed_concord ~= "0");
% 
%     ContTable_2x2 = [n0_disc,  nNon0_disc;
%                      n0_conc, nNon0_conc];
% 
%     if sum(ContTable_2x2(1,:))>0 && sum(ContTable_2x2(2,:))>0
%         % Then we can do a normal 2×2 fisher
%         [~, pVal] = fishertest(ContTable_2x2);
%         fisherP_vals(i) = pVal;
%     else
%         fisherP_vals(i) = NaN;
%     end
% 
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('Contingency table 2×2 (Discord vs. Conc) x (AED=0 vs. AED=1):');
%     disp(ContTable_2x2);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
% end
% 
% % ----------- Create a Stacked Bar for DISCORDANT Only -------------
% barDataAED = [n0, n1, n2, n3plus];   % [nROIs x 4]
% 
% figure('Color','w','Name','Discordant By AED count','Position',[100 300 700 400]);
% hBar = bar(barDataAED, 'stacked');
% 
% % Optional color scheme
% colorSpec = [ ...
%     0.8 0.2 0.2;   % '0'
%     0.5 0.5 0.5;   % '1'
%     0.2 0.6 0.8;   % '2'
%     0.4 0.2 0.6];  % '3plus'
% for c = 1:size(barDataAED,2)
%     hBar(c).FaceColor = colorSpec(c,:);
% end
% 
% set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize', 10);
% xlabel('ROI');
% ylabel('# Discordant Patients');
% title('Discordant Subjects by AED Count');
% legend({'0','1','2','3plus'}, 'Location','bestoutside');
% box off;
% 
% % Print final p-values
% for i = 1:nROIs
%     fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
% end


% allROI = bestResultsTable.ROI;       % cell array of ROI names
% nROIs  = height(bestResultsTable);   % number of rows (ROIs)
% 
% allSubs = 1:72;  % total # subjects (adjust if needed)
% 
% % Pre-allocate arrays to store counts for each AED bin (DISCORDANT only)
% n0      = zeros(nROIs,1);
% n1      = zeros(nROIs,1);
% n2      = zeros(nROIs,1);
% n3plus  = zeros(nROIs,1);
% 
% % Store multi-way Fisher p-values:
% fisherP_vals = nan(nROIs,1);
% 
% for i = 1:nROIs
% %     discordSubs = bestResultsTable.Best_Discord_Subs{i};
%     discordSubs = bestResultsTable.Gross_Discord_Subs{i};
% 
%     concordSubs = setdiff(allSubs, discordSubs);
%     
%     % Subset the AED categories
%     aed_discord = T1_epil_measures_upted.AEDcat(discordSubs);
%     aed_concord = T1_epil_measures_upted.AEDcat(concordSubs);
%     
%     % Count '0','1','2','3plus' for DISCORDANT
%     n0(i)     = sum(aed_discord == '0');
%     n1(i)     = sum(aed_discord == '1');
%     n2(i)     = sum(aed_discord == '2');
%     n3plus(i) = sum(aed_discord == '3plus');
%     
%     % 2×4 Table:
%     n0_disc = sum(aed_discord == '0');
%     n1_disc = sum(aed_discord == '1');
%     n2_disc = sum(aed_discord == '2');
%     n3_disc = sum(aed_discord == '3plus');
%     
%     n0_conc = sum(aed_concord == '0');
%     n1_conc = sum(aed_concord == '1');
%     n2_conc = sum(aed_concord == '2');
%     n3_conc = sum(aed_concord == '3plus');
%     
%     ContTable = [n0_disc, n1_disc, n2_disc, n3_disc;
%         n0_conc, n1_conc, n2_conc, n3_conc];
%     
%     % If each row has at least one total subject, we try fishertest
%     if sum(ContTable(1,:))>0 && sum(ContTable(2,:))>0
%         %         [~, pVal] = fishertest(ContTable, 'Tail','both');
%         n0_disc   = sum(aed_discord == '0');
%         nNon0_disc= sum(aed_discord ~= '0');
%         n0_conc   = sum(aed_concord == '0');
%         nNon0_conc= sum(aed_concord ~= '0');
%         
%         ContTable_2x2 = [n0_disc, nNon0_disc; n0_conc, nNon0_conc];
%         [~, pVal] = fishertest(ContTable_2x2);
%         fisherP_vals(i) = pVal;
%     else
%         fisherP_vals(i) = NaN;
%     end
%     
%     % Print out each ROI's table & p-value
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('2×4 Fisher test (Discordant/Concordant × {0,1,2,3plus}):');
%     disp(ContTable);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
% end
% 
% % Create a stacked bar for the DISCORDANT group
% barDataAED = [n0, n1, n2, n3plus];
% figure('Color','w','Name','Discordant By AED count','Position',[100 300 700 400]);
% hBar = bar(barDataAED, 'stacked');
% 
% colorSpec = [ ...
%     0.8 0.2 0.2;   % '0'
%     0.5 0.5 0.5;   % '1'
%     0.2 0.6 0.8;   % '2'
%     0.4 0.2 0.6];  % '3plus'
% for c = 1:size(barDataAED,2)
%     hBar(c).FaceColor = colorSpec(c,:);
% end
% 
% set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize', 10);
% xlabel('ROI');
% ylabel('# Discordant Patients');
% title('Discordant Subjects by AED Count');
% legend({'0','1','2','3plus'}, 'Location','bestoutside');
% box off;
% 
% % Print final p-values
% for i = 1:nROIs
%     fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
% end
