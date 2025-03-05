%% run_plot_CPFreq_simplified.m
%
% We merge 'cp_freq_cat' into 2 bins:
%   '0'   vs.   '>=1'  (where '>=1' lumps '1to5','6to10','>10')
%
% Then we do a 2×2 Fisher test: (Discordant vs Concordant) × (CP=0 vs CP=1)
% and plot a 2-column stacked bar for the Discordant group.

%% A) Basic Setup
allROI = bestResultsTable.ROI;      % cell array of ROI names
nROIs  = height(bestResultsTable);  % number of rows
allSubs = 1:72;                     % total # subjects (adjust if needed)

% Decide if you use "Best_Discord_Subs" or "Gross_Discord_Subs"
discordColumn = 'Best_Discord_Subs';
% or 'Gross_Discord_Subs' if you want "gross mismatch"

%% B) We'll store two columns: col1='0', col2='=1'
n_zero = zeros(nROIs,1);
n_plus = zeros(nROIs,1);

% We'll store p-values from the 2×2 Fisher
fisherP_vals = nan(nROIs,1);

for i = 1:nROIs
    %% 1) Identify Discordant vs. Concordant
    discordSubs = bestResultsTable.(discordColumn){i};
    concordSubs = setdiff(allSubs, discordSubs);

    %% 2) Subset T1_epil_measures_upted.cp_freq_cat
    cp_discord = T1_epil_measures_upted.cp_freq_cat(discordSubs);
    cp_concord = T1_epil_measures_upted.cp_freq_cat(concordSubs);

    %% 3) For DISCORDANT: count how many are '0' vs '=1'
    n_zero(i) = sum(cp_discord == "0");
    n_plus(i) = sum(cp_discord ~= "0");  % lumps '1to5','6to10','>10'

    %% 4) 2×2 Fisher: (Discord vs Concord) × (CP=0 vs CP=1)
    zero_disc   = sum(cp_discord == "0");
    plus_disc   = sum(cp_discord ~= "0");
    zero_conc   = sum(cp_concord == "0");
    plus_conc   = sum(cp_concord ~= "0");

    ContTable_2x2 = [zero_disc, plus_disc; zero_conc, plus_conc];

    if sum(ContTable_2x2(1,:))>0 && sum(ContTable_2x2(2,:))>0
        [~, pVal] = fishertest(ContTable_2x2);
        fisherP_vals(i) = pVal;
    else
        fisherP_vals(i) = NaN;
    end

%     % Display
%     fprintf('\nROI #%d = %s\n', i, allROI{i});
%     disp('(Discord vs Concord) × (CP=0 vs CP>=1) :');
%     disp(ContTable_2x2);
%     fprintf('p-value = %.4g\n', fisherP_vals(i));
end

%% C) Create the 2-column stacked bar for DISCORDANT
barData = [n_zero, n_plus];  % Nx2

figure('Color','w','Name','Discordant By CP_Freq (0 vs >=1)','Position',[100 300 700 400]);
hBar = bar(barData, 'stacked');

% Color scheme for 2 columns
colorSpec = [0.8 0.2 0.2; 0.2 0.6 0.8];  % red, teal
for c = 1:size(barData,2)
    hBar(c).FaceColor = colorSpec(c,:);
end

set(gca,'XTick',1:nROIs,'XTickLabel',allROI,'FontSize',9);
xlabel('ROI');
ylabel('# Discordant Patients');
title('Discordant Subjects by CP Frequency (0 vs >=1)');
legend({'0','>=1'}, 'Location','bestoutside');
box off;

%% D) Print final p-values
for i = 1:nROIs
    fprintf('ROI: %s | p=%.4g\n', allROI{i}, fisherP_vals(i));
end


% %% run_plot_CPFreq.m
% 
% allROI = bestResultsTable.ROI;
% nROIs  = height(bestResultsTable);
% 
% CPcats = categories(T1_epil_measures_upted.cp_freq_cat);  % e.g. {'0','1to5','6to10','>10'}
% nCats  = length(CPcats);
% 
% CP_counts = zeros(nROIs, nCats);
% 
% for i = 1:nROIs
% %     discordSubs = bestResultsTable.Best_Discord_Subs{i};
%     discordSubs = bestResultsTable.Gross_Discord_Subs{i};
%     cpVals      = T1_epil_measures_upted.cp_freq_cat(discordSubs);
% 
%     for c = 1:nCats
%         CP_counts(i,c) = sum(cpVals == CPcats{c});
%     end
% end
% 
% figure('Color','w','Name','Discordant By CP_Freq','Position',[100 300 700 400]);
% hBar = bar(CP_counts, 'stacked');
% 
% colorMap = [ ...
%     0.8 0.2 0.2;  ... red
%     0.5 0.5 0.5;  ... gray
%     0.2 0.6 0.8;  ... teal
%     0.4 0.2 0.6;  ... more if you have more categories
% ];
% for c = 1:min(nCats,size(colorMap,1))
%     hBar(c).FaceColor = colorMap(c,:);
% end
% 
% set(gca,'XTick',1:nROIs,'XTickLabel',allROI,'FontSize',9);
% xlabel('ROI');
% ylabel('# Discordant Patients');
% title('Discordant Subjects by CP Frequency');
% legend(CPcats, 'Location','bestoutside');
% box off;


