%% run_plot_CPfreq.m
allROI = bestResultsTable.ROI;
nROIs  = height(bestResultsTable);

n0     = zeros(nROIs,1);
n1to5  = zeros(nROIs,1);
n6to10 = zeros(nROIs,1);
n11p   = zeros(nROIs,1);

for i = 1:nROIs
    discordSubs = bestResultsTable.Best_Discord_Subs{i};
    cp_discord  = T1_epil_measures.cp_freq_cat(discordSubs);

    n0(i)     = sum(cp_discord == '0');
    n1to5(i)  = sum(cp_discord == '1to5');
    n6to10(i) = sum(cp_discord == '6to10');
    n11p(i)   = sum(cp_discord == '11plus');
end

figure('Color','w','Name','Discordant By CP_Freq','Position',[100 300 700 400]);
barData = [n0, n1to5, n6to10, n11p];
hBar = bar(barData, 'stacked');

colors = [0.8 0.2 0.2;  ... red
          0.5 0.5 0.5;  ... gray
          0.2 0.6 0.8;  ... teal
          0.6 0.4 0.8]; ... purple
for c = 1:4
    hBar(c).FaceColor = colors(c,:);
end

set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize', 9);
xlabel('ROI');
ylabel('# Discordant Patients');
title('Discordant Subjects by CP-Freq Bins');
legend({'0','1-5','6-10','11+'}, 'Location','bestoutside');
box off;
