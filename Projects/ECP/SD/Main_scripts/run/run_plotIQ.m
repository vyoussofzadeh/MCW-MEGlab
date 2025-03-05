allROI = bestResultsTable.ROI;       % cell array of ROI names
nROIs  = height(bestResultsTable);   % number of rows (ROIs)

allSubs = 1:72;  % total subjects (adjust if needed)

% Pre-allocate for the stacked bar of DISCORDANT only
nAboveIQ  = zeros(nROIs,1);
nAvgIQ    = zeros(nROIs,1);
nBelowIQ  = zeros(nROIs,1);

% (Optional) store Fisher p-values and odds ratios
fisherP_vals  = nan(nROIs,1);
oddsRatio_vals= nan(nROIs,1);

for i = 1:nROIs
    % 1) Identify Discordant vs. Concordant for this ROI
%     discordSubs = bestResultsTable.Best_Discord_Subs{i};
    discordSubs = bestResultsTable.Gross_Discord_Subs{i};

    concordSubs = setdiff(allSubs, discordSubs);

    % 2) Extract IQ categories for DISCORDANT subset
    IQcat_discord  = T1_epil_measures_upted.IQcat(discordSubs);

    % A) Count how many are Above, Average, Below for discordant 
    nAboveIQ(i) = sum(IQcat_discord == "Above");
    nAvgIQ(i)   = sum(IQcat_discord == "Average");
    nBelowIQ(i) = sum(IQcat_discord == "Below");

    % ------ Fishers Exact (2×2) ignoring 'Average' ------
    % We'll only compare 'Above' vs. 'Below' for both groups.

    % (a) Discordant counts
    nAboveDiscord = sum(IQcat_discord == "Above");
    nBelowDiscord = sum(IQcat_discord == "Below");
    nAverageDiscord = sum(IQcat_discord == "Average");
    
    % (b) Concordant counts
    IQcat_concord = T1_epil_measures_upted.IQcat(concordSubs);
    nAboveConcord = sum(IQcat_concord == "Above");
    nBelowConcord = sum(IQcat_concord == "Below");
    nAverageConcord = sum(IQcat_concord == "Average");


    % (c) Build 2×2: rows = (Discordant,Concordant), cols = (Above, Below)
    ContTable = [nAverageDiscord,  nBelowDiscord;
                 nAverageConcord, nBelowConcord];

    % If a row is all zeros or something, fishertest might fail
    % We'll skip the test in that case:
    if sum(ContTable(1,:))>0 && sum(ContTable(2,:))>0
        [~, pVal, stats] = fishertest(ContTable);
        fisherP_vals(i)   = pVal;
        if isfield(stats,'OddsRatio')
            oddsRatio_vals(i) = stats.OddsRatio;
        end
    else
        fisherP_vals(i)   = NaN;
        oddsRatio_vals(i) = NaN;
    end
    
%     % Print results for each ROI
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('Fisher Exact Test (IQ Above vs. Below) by Concordance vs. Discordance:');
%     disp(ContTable);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
%     if ~isnan(oddsRatio_vals(i))
%         fprintf('OddsRatio = %.3f\n', oddsRatio_vals(i));
%     end
end

% ------------ PLOT: Stacked Bar of DISCORDANT ONLY -------------
%   Each bar has [Above, Average, Below] for that ROI.
figure('Color','w','Name','Discordant By IQ Category','Position',[100 300 700 400]);
barDataIQ = [nAboveIQ, nAvgIQ, nBelowIQ];
hBar2 = bar(barDataIQ, 'stacked');

% Color each segment
hBar2(1).FaceColor = [0.8 0.2 0.2];  % red (Above)
hBar2(2).FaceColor = [0.5 0.5 0.5];  % gray (Average)
hBar2(3).FaceColor = [0.2 0.6 0.8];  % teal (Below)

set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize', 10);
xlabel('ROI');
ylabel('# Discordant Patients');
title('Discordant Subjects by IQ Category');
legend({'Above','Average','Below'}, 'Location','bestoutside');
box off;

% ---------(Optional) Check p-values across ROIs-----------
for i = 1:nROIs
    fprintf('ROI: %s | Fisher p=%.4g | OR=%.3f\n', ...
        allROI{i}, fisherP_vals(i), oddsRatio_vals(i));
end


% nAboveIQ = zeros(nROIs,1);
% nAvgIQ   = zeros(nROIs,1);
% nBelowIQ = zeros(nROIs,1);
% 
% for i = 1:nROIs
%     discordSubs = bestResultsTable.Best_Discord_Subs{i};
%     IQcat_discord = T1_epil_measures_upted.IQcat(discordSubs);
% 
%     nAboveIQ(i)  = sum(IQcat_discord == 'Above');
%     nAvgIQ(i)    = sum(IQcat_discord == 'Average');
%     nBelowIQ(i)  = sum(IQcat_discord == 'Below');
% end
% 
% figure('Color','w','Name','Discordant By IQ Category','Position',[100 300 700 400]);
% barDataIQ = [nAboveIQ, nAvgIQ, nBelowIQ];
% hBar2 = bar(barDataIQ, 'stacked');
% 
% % Reuse the color scheme or pick new ones:
% hBar2(1).FaceColor = [0.8 0.2 0.2];  % red
% hBar2(2).FaceColor = [0.5 0.5 0.5];  % gray
% hBar2(3).FaceColor = [0.2 0.6 0.8];  % teal
% 
% set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize', 10);
% xlabel('ROI');
% ylabel('# Discordant Patients');
% title('Discordant Subjects by IQ');
% legend({'Above IQ','Average IQ','Below IQ'}, 'Location','bestoutside');
% box off;
