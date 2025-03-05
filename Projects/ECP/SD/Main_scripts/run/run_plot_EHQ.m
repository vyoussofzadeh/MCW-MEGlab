% Pre-allocate for the stacked bar (DISCORDANT) by EHQ categories
nLeftEHQ  = zeros(nROIs,1);
nAmbiEHQ  = zeros(nROIs,1);
nRightEHQ = zeros(nROIs,1);

% (Optional) store Fisher's p-values & odds ratios
fisherP_vals  = nan(nROIs,1);
oddsRatio_vals= nan(nROIs,1);

allSubs = 1:72;  % total # subjects (adjust as needed)

for i = 1:nROIs
    % 1) Identify Discordant vs. Concordant for this ROI
%     discordSubs = bestResultsTable.Best_Discord_Subs{i};
    discordSubs = bestResultsTable.Gross_Discord_Subs{i};

    concordSubs = setdiff(allSubs, discordSubs);

    % 2) Extract EHQ categories for the discordant subset
    ehqCat_discord = T1_epil_measures_upted.EHQcat(discordSubs);

    % A) Count how many are Left, Ambi, Right among discordant
    nLeftEHQ(i)  = sum(ehqCat_discord == "Left");
    nAmbiEHQ(i)  = sum(ehqCat_discord == "Ambi");
    nRightEHQ(i) = sum(ehqCat_discord == "Right");

    % ----- 2×2 Fisher ignoring 'Ambi' -----
    nLeftDiscord  = sum(ehqCat_discord == "Left");
    nRightDiscord = sum(ehqCat_discord == "Right");

    ehqCat_concord = T1_epil_measures_upted.EHQcat(concordSubs);
    nLeftConcord   = sum(ehqCat_concord == "Left");
    nRightConcord  = sum(ehqCat_concord == "Right");

    ContTable = [nLeftDiscord,  nRightDiscord; 
                 nLeftConcord, nRightConcord];

    % Ensure each row > 0
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
    
%     % Print results
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('Fisher Exact Test: Left vs Right by Concordance vs. Discordance');
%     disp(ContTable);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
%     if ~isnan(oddsRatio_vals(i))
%         fprintf('OddsRatio = %.3f\n', oddsRatio_vals(i));
%     end
end

%% Plot stacked bar of DISCORDANT subjects (Left, Ambi, Right)
figure('Color','w','Name','Discordant By EHQ Category','Position',[100 300 700 400]);
barDataEHQ = [nLeftEHQ, nAmbiEHQ, nRightEHQ];
hBar3 = bar(barDataEHQ, 'stacked');

% Color scheme
hBar3(1).FaceColor = [0.8 0.2 0.2];  % red (Left)
hBar3(2).FaceColor = [0.5 0.5 0.5];  % gray (Ambi)
hBar3(3).FaceColor = [0.2 0.6 0.8];  % teal (Right)

set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize', 10);
xlabel('ROI');
ylabel('# Discordant Subjects');
title('Discordant Subjects by Handedness (EHQ)');
legend({'Left-Handed','Ambidextrous','Right-Handed'}, 'Location','bestoutside');
box off;

%% Inspect the Fisher p-values across all ROIs
for i = 1:nROIs
    fprintf('ROI: %s | p=%.4g | OR=%.3f\n', allROI{i}, fisherP_vals(i), oddsRatio_vals(i));
end


% nLeftEHQ  = zeros(nROIs,1);
% nAmbiEHQ  = zeros(nROIs,1);
% nRightEHQ = zeros(nROIs,1);
% 
% for i = 1:nROIs
%     discordSubs = bestResultsTable.Best_Discord_Subs{i};
%     ehqCat_discord = T1_epil_measures_upted.EHQcat(discordSubs);
% 
%     nLeftEHQ(i)  = sum(ehqCat_discord == 'Left');
%     nAmbiEHQ(i)  = sum(ehqCat_discord == 'Ambi');
%     nRightEHQ(i) = sum(ehqCat_discord == 'Right');
% end
% 
% figure('Color','w','Name','Discordant By EHQ Category','Position',[100 300 700 400]);
% barDataEHQ = [nLeftEHQ, nAmbiEHQ, nRightEHQ];
% hBar3 = bar(barDataEHQ, 'stacked');
% 
% hBar3(1).FaceColor = [0.8 0.2 0.2];  % red
% hBar3(2).FaceColor = [0.5 0.5 0.5];  % gray
% hBar3(3).FaceColor = [0.2 0.6 0.8];  % teal
% 
% set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize', 10);
% xlabel('ROI');
% ylabel('# Discordant Subjects');
% title('Discordant Subjects by Handedness (EHQ)');
% legend({'Left-Handed','Ambidextrous','Right-Handed'}, 'Location','bestoutside');
% box off;
% 
