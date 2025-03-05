%% run_plot_LTGTC.m
%
% This script takes:
%   - bestResultsTable: with columns
%       * Best_Discord_Subs (or Gross_Discord_Subs) for each ROI
%   - T1_epil_measures: with a column 'LTGTCcat' (categorical)
%   - N subjects total (e.g. 72)
%
% We merge all categories like {"0","1","2-3","4-5","6-10","11-20"}
% into a single bin "<=20",
% and everything else {"21-50","51-100",">100"} into ">20".
%
% Then we do a 2x2 Fisher test: (Discordant vs. Concordant) x (<=20 vs. >20)
% Finally, we plot a 2-bin stacked bar for the Discordant group only.

%% (1) Basic Setup
allROI = bestResultsTable.ROI;  % cell array of ROI names
nROIs  = height(bestResultsTable);

allSubs = 1:72;  % total # of subjects (adjust if needed)

% Decide which "discordant" we use: 
%   bestResultsTable.Best_Discord_Subs{i}  or  bestResultsTable.Gross_Discord_Subs{i}
discordColumn = 'Best_Discord_Subs';
% If you want "gross" mismatch, do: discordColumn = 'Gross_Discord_Subs';

%% (2) Define bigger bins: "=20" vs. ">20"
% We'll create two sets of categories (both as categorical arrays) so we can do ismember( , ) properly.
le20Cats = categorical(["0","1","2-3","4-5","6-10","11-20"]);    % =20
gt20Cats = categorical(["21-50","51-100",">100"]);              % >20

% We'll create inline functions to check membership in those sets
isLE20 = @(catArray) ismember(catArray, le20Cats);
isGT20 = @(catArray) ismember(catArray, gt20Cats);

%% Pre-allocate counters for the stacked bar (2 columns => "=20", ">20") for DISCORDANT only
n_le20 = zeros(nROIs,1);
n_gt20 = zeros(nROIs,1);

% We'll store p-values from the Fisher test in an array
fisherP_vals = nan(nROIs,1);

for i = 1:nROIs
    % (A) Get the "discordant" subset for this ROI
    discordSubs = bestResultsTable.(discordColumn){i};  % e.g. Best_Discord_Subs{i}
    % The rest are "concordant"
    concordSubs = setdiff(allSubs, discordSubs);

    % (B) Extract the LTGTCcat for these subsets
    LTdiscord = T1_epil_measures_upted.LTGTCcat(discordSubs);
    LTconcord= T1_epil_measures_upted.LTGTCcat(concordSubs);

    % (C) For the DISCORDANT group, count how many are in the "=20" bin vs. ">20"
    n_le20(i) = sum(isLE20(LTdiscord));
    n_gt20(i) = sum(isGT20(LTdiscord));
    % If there are categories like 'Missing' or 'Unknown', they won't count in either bin.

    % -------- 2×2 Fisher: (Discord vs. Concord) × (=20 vs. >20) --------
    le20_disc = sum(isLE20(LTdiscord));
    gt20_disc = sum(isGT20(LTdiscord));
    le20_conc = sum(isLE20(LTconcord));
    gt20_conc = sum(isGT20(LTconcord));

    ContTable_2x2 = [le20_disc,  gt20_disc;
                     le20_conc, gt20_conc];

    % (D) If each row > 0, attempt fishertest
    if sum(ContTable_2x2(1,:))>0 && sum(ContTable_2x2(2,:))>0
        [~, pVal] = fishertest(ContTable_2x2);
        fisherP_vals(i) = pVal;
    else
        fisherP_vals(i) = NaN;
    end

%     % (E) Display results
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('ContTable (Discord vs Concord) x (<=20 vs >20):');
%     disp(ContTable_2x2);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
end

%% (3) Plot the 2-column stacked bar for DISCORDANT only
barData = [n_le20, n_gt20];  % [nROIs x 2]

figure('Color','w','Name','Discordant By LTGTC(BiggerBins)','Position',[100 300 700 400]);
hBar = bar(barData, 'stacked');

% Optional color scheme: 2 columns
colorSpec = [0.2 0.6 0.8;  0.8 0.2 0.2];  % e.g. teal, red
for c = 1:size(barData,2)
    hBar(c).FaceColor = colorSpec(c,:);
end

set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize',10);
xlabel('ROI');
ylabel('# Discordant Patients');
title('Discordant Subjects by LTGTC (<=20 vs >20)');

legend({'<=20','>20'}, 'Location','bestoutside');
box off;

%% (4) Print final p-values for convenience
for i = 1:nROIs
    fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
end



% %% run_plot_LTGTC.m
% % Plots how many DISCORDANT subjects are in each LTGTCcat 
% % category (c_0 => '0', c_4-5 => '4-5', etc.) across ROIs, 
% % and does a 2×n Fisher test (Discordant vs Concordant) × (LTGTC categories).
% 
% % (1) Basic setup
% allROI = bestResultsTable.ROI;
% nROIs  = height(bestResultsTable);
% 
% % Make sure you know the total # of subjects
% allSubs = 1:72;  % or however many subjects you have
% 
% % (2) We get the categories from T1_epil_measures_upted.LTGTCcat
% %    e.g. {'0','1','2-3','4-5','6-10','21-50',...}
% LTcats = categories(T1_epil_measures_upted.LTGTCcat);  
% nCats  = length(LTcats);
% 
% % (3) Build a [nROIs x nCats] matrix for the DISCORDANT group
% LT_counts = zeros(nROIs, nCats);
% 
% % (Optional) store p-values from the multi-way Fisher test for each ROI
% fisherP_vals = nan(nROIs,1);
% 
% for i = 1:nROIs
%     % -- 3A) Identify discordant vs. concordant for this ROI
% %     discordSubs  = bestResultsTable.Best_Discord_Subs{i};
%     discordSubs = bestResultsTable.Gross_Discord_Subs{i};
%     concordSubs  = setdiff(allSubs, discordSubs);
% 
%     % -- 3B) LTGTC categories for discordant
%     theseCatsDiscord  = T1_epil_measures_upted.LTGTCcat(discordSubs);
% 
%     % Count how many subjects in each category of LTcats
%     for c = 1:nCats
%         LT_counts(i, c) = sum(theseCatsDiscord == LTcats{c});
%     end
% 
%     % ---- 4) Fishers exact (2×n) ignoring row/column if sum=0
%     if ~isempty(discordSubs) && ~isempty(concordSubs)
%         % Build the 2×n table
%         %  Row1 = Discordant distribution among all categories
%         %  Row2 = Concordant distribution among all categories
%         ContTable = zeros(2, nCats);
% 
%         for c = 1:nCats
%             ContTable(1,c) = sum(theseCatsDiscord == LTcats{c});
%             ContTable(2,c) = sum(T1_epil_measures_upted.LTGTCcat(concordSubs) == LTcats{c});
%         end
% 
%         % If each row has at least 1 subject, attempt fishertest
%         if sum(ContTable(1,:))>0 && sum(ContTable(2,:))>0
%             try
%                 % 2×n Fisher
%                 [~, pVal] = fishertest(ContTable,'Tail','both');
%                 fisherP_vals(i) = pVal;
%             catch ME
%                 % If fishertest errors (older MATLAB), set pVal=NaN
%                 warning('Fisher test error in ROI %s: %s', allROI{i}, ME.message);
%                 fisherP_vals(i) = NaN;
%             end
%         end
%     end
% 
%     % Print results
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('2×n Fisher test (Discordant/Concordant × LTGTC categories):');
%     disp('ContTable =');
%     disp(ContTable);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
% end
% 
% %% (5) Plot stacked bar for DISCORDANT only
% figure('Color','w','Name','Discordant By LTGTC','Position',[100 300 700 400]);
% hBar = bar(LT_counts, 'stacked');
% 
% % (6) Optional color scheme (with up to 10 rows).
% colorMap = [ ...
%     0.8 0.2 0.2;   % '0'
%     0.5 0.5 0.5;   % '1'
%     0.2 0.6 0.8;   % '2-3'
%     0.6 0.4 0.8;   % '4-5'
%     0.9 0.6 0.2;   % '6-10'
%     0.2 0.8 0.2;   % '11-20'
%     0.8 0.2 0.4;   % '21-50'
%     0.3 0.5 0.3;   % '51-100'
%     0.4 0.2 0.6;   % '>100'
%     0.7 0.7 0.7];  % 'Missing'
% 
% for c = 1:min(nCats, size(colorMap,1))
%     hBar(c).FaceColor = colorMap(c,:);
% end
% 
% % (7) Labeling
% set(gca, 'XTick',1:nROIs, 'XTickLabel', allROI, 'FontSize',10);
% xlabel('ROI');
% ylabel('# Discordant Patients');
% title('Discordant Subjects by LTGTC (Multi-level)');
% 
% legend(LTcats, 'Location','bestoutside');
% box off;
% 
% %% Print final p-values
% for i = 1:nROIs
%     fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
% end
