% Initialize a table to store the best results
bestResultsTable_fixed = table();

uniqueROIs = unique(summaryTable.ROI);


% Loop through each ROI and select the best LI method based on highest Concordance
for i = 1:length(uniqueROIs)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = summaryTableFixed(strcmp(summaryTableFixed.ROI, roi), :);
    
    % Initialize variables to keep track of the best results
    bestCorrelation = -inf;
    bestConcordance = -inf;
    bestMethod = '';
    bestTimeInterval = [];
    bestDiscordSubs = [];
    
    % Loop through each method and find the best one based on Concordance
    for j = 1:height(roiData)
        methodData = roiData(j, :);
        
        if methodData.Concordance > bestConcordance
            bestConcordance = methodData.Concordance;
            bestCorrelation = methodData.Correlation;
            bestMethod = methodData.LI_Method;
            Corr_P_value = methodData.Corr_P_value;
            kappa = methodData.kappa;
            bestTimeInterval = methodData.mean_Optimal_Time;
            bestDiscordSubs = methodData.discord_Subs;
            bestDiscordSubs_megli = methodData.MEG_LI;
            bestDiscordSubs_optmegli = methodData.optimalMEG_LI;
            bestDiscordSubs_fmrili = methodData.fMRI_LI;
            bestDiscordSubs_rSNR_L = methodData.rSNR_left;
            bestDiscordSubs_rSNR_R = methodData.rSNR_right;            
        end
    end
    
    % Store the best results in the table
    newRow = {bestMethod, roi, bestCorrelation, Corr_P_value, bestConcordance, kappa, bestTimeInterval, bestDiscordSubs, bestDiscordSubs_megli, bestDiscordSubs_optmegli, bestDiscordSubs_fmrili, bestDiscordSubs_rSNR_L, bestDiscordSubs_rSNR_R};
    bestResultsTable_fixed = [bestResultsTable_fixed; newRow];
end

% Set column names for the best results table
bestResultsTable_fixed.Properties.VariableNames = {'Best_LI_Method', 'ROI', 'Best_Correlation', 'Corr_pval', 'Best_Concordance', 'kappa', 'Best_Time_Interval', 'Best_Discord_Subs', 'MEG_LI', 'optMEG_LI', 'fMRI_LI', 'rSNR_L', 'rSNR_R'};

bestResultsTable_save = bestResultsTable_fixed;
bestResultsTable_save = bestResultsTable_save(:,1:end-3);



% Format the numeric data to two decimal places before saving to CSV
bestResultsTable_save.Best_Correlation = round(bestResultsTable_save.Best_Correlation, 2);
bestResultsTable_save.Best_Concordance = round(bestResultsTable_save.Best_Concordance, 2);
bestResultsTable_save.Best_Time_Interval = round(bestResultsTable_save.Best_Time_Interval, 2);

% Save best results table
writetable(bestResultsTable_save, fullfile(save_dir,'Best_LI_Methods_Summary.csv'));


% Save best results table as text file
fid = fopen(fullfile(save_dir,'Best_LI_Methods_Summary.txt'), 'wt');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', bestResultsTable_fixed.Properties.VariableNames{:});

for i = 1:height(bestResultsTable_fixed)
    fprintf(fid, '%s\t%s\t%.2f\t%.2f\t%.2f\n', bestResultsTable_fixed.Best_LI_Method{i}, bestResultsTable_fixed.ROI{i}, ...
        bestResultsTable_fixed.Best_Correlation(i), bestResultsTable_fixed.Best_Concordance(i), bestResultsTable_fixed.Best_Time_Interval(i));
end

fclose(fid);

disp('Best LI methods for discordant analyses.');
disp('bestResultsTable_fixed:')
disp(bestResultsTable_fixed);