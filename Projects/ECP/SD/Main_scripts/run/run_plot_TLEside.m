allROI = bestResultsTable.ROI;       % cell array of ROI names
nROIs  = height(bestResultsTable);   % number of rows (ROIs)

% Pre-allocate counters for TLE side (only for discordant bars)
nLeftAll_discord  = zeros(nROIs,1);
nRightAll_discord = zeros(nROIs,1);
nBilatAll_discord = zeros(nROIs,1);

% (Optional) store the Fisher p-value for each ROI
fisherP_vals  = zeros(nROIs,1);
oddsRatio_vals= zeros(nROIs,1);

for i = 1:nROIs
    % (A) Grab discordant subject indices for this ROI
%     discordSubs = bestResultsTable.Best_Discord_Subs{i};
    discordSubs = bestResultsTable.Gross_Discord_Subs{i};

    % Suppose you have 72 total subjects
    allSubs = 1:72;  
    concordSubs = setdiff(allSubs, discordSubs);
    
    % (B) Get TLE side for these subsets
    TLEside_discord = T1_epil_measures_upted.TLEside(discordSubs);  % 'Left','Right','Bilateral'
    TLEside_concord = T1_epil_measures_upted.TLEside(concordSubs);
    
    % (C) Count how many have Left, Right, Bilateral *for the DISCORDANT group*
    nLeftAll_discord(i)   = sum(TLEside_discord == 'Left');
    nRightAll_discord(i)  = sum(TLEside_discord == 'Right');
    nBilatAll_discord(i)  = sum(TLEside_discord == 'Bilateral');
    
    % If you want the total # of Concordant left, right, bilateral as well,
    % you can store them similarly in nLeftAll_concord, etc.

    % ----- 2x2 Fishers exact ignoring Bilateral -----
    % Count only left & right for Discordant
    nLeftDiscord  = sum(TLEside_discord == 'Left');
    nRightDiscord = sum(TLEside_discord == 'Right');
    
    % Count only left & right for Concordant
    nLeftConcord  = sum(TLEside_concord == 'Left');
    nRightConcord = sum(TLEside_concord == 'Right');
    
    % Build the 2x2 table: rows = (Discordant, Concordant), cols = (Left, Right)
    ContTable = [nLeftDiscord,  nRightDiscord;
                 nLeftConcord, nRightConcord];
    
    % If either row sums to 0, fishertest can fail. Check or skip if needed:
    if all(sum(ContTable,2) > 0)
        [~, fisherP, stats] = fishertest(ContTable);
    else
        % If there's a row of zeros, the test won't run properly
        fisherP = NaN;
        stats.OddsRatio = NaN;
    end
    fisherP_vals(i)   = fisherP;
    oddsRatio_vals(i) = stats.OddsRatio;
    
    % Print or store results
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('Fisher Exact Test: TLE=Left vs Right by Concordance/Discordance');
%     disp(ContTable);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
    
%     if ~isnan(oddsRatio_vals(i))
%         fprintf('Odds Ratio = %.3f\n', oddsRatio_vals(i));
%         if isfield(stats,'ConfidenceInterval')
%             disp('95% CI:'); disp(stats.ConfidenceInterval);
%         end
%     end
end

% ---------- Now create a stacked bar of # DISCORDANT subjects by TLE side ----------
% each row = ROI, columns = [Left, Bilateral, Right]
figure('Color','w','Name','Discordant By TLE Side','Position',[100 300 700 400]);
barData = [nLeftAll_discord, nBilatAll_discord, nRightAll_discord];
hBar = bar(barData, 'stacked');

% Assign colors
hBar(1).FaceColor = [0.8 0.2 0.2];  % red
hBar(2).FaceColor = [0.5 0.5 0.5];  % gray
hBar(3).FaceColor = [0.2 0.6 0.8];  % teal

set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize', 10);
xlabel('ROI');
ylabel('# Discordant Patients');
title('Discordant Subjects by TLE Side');
legend({'Left','Bilateral','Right'}, 'Location','bestoutside');
box off;

% ---------- (Optional) Inspect the Fisher p-values for each ROI ----------
for i = 1:nROIs
    fprintf('ROI: %s | Fisher p=%.4f | OR=%.3f\n', ...
             allROI{i}, fisherP_vals(i), oddsRatio_vals(i));
end
