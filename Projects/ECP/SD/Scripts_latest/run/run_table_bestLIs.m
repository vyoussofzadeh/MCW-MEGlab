% Initialize a table to store the best results
bestResultsTable = table();

% Loop through each ROI and select the best LI method based on highest Concordance
for i = 1:length(uniqueROIs)
    roi = uniqueROIs{i};
    
    % Extract data for the current ROI
    roiData = summaryTableDynamic(strcmp(summaryTableDynamic.ROI, roi), :);
    
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
            bestTimeInterval = methodData.mean_Optimal_Time;
            bestDiscordSubs = methodData.discord_Subs;
            bestDiscordSubs_megli = methodData.MEG_LI;
            bestDiscordSubs_fmrili = methodData.fMRI_LI;
        end
    end
    
    % Store the best results in the table
    newRow = {bestMethod, roi, bestCorrelation, bestConcordance, bestTimeInterval, bestDiscordSubs, bestDiscordSubs_megli, bestDiscordSubs_fmrili};
    bestResultsTable = [bestResultsTable; newRow];
end

% Set column names for the best results table
bestResultsTable.Properties.VariableNames = {'Best_LI_Method', 'ROI', 'Best_Correlation', 'Best_Concordance', 'Best_Time_Interval', 'Best_Discord_Subs', 'MEG_LI', 'fMRI_LI'};

% Save best results table
writetable(bestResultsTable, 'Best_LI_Methods_Summary.csv');

% Save best results table as text file
fid = fopen('Best_LI_Methods_Summary.txt', 'wt');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', bestResultsTable.Properties.VariableNames{:});

for i = 1:height(bestResultsTable)
    fprintf(fid, '%s\t%s\t%f\t%f\t%f\n', bestResultsTable.Best_LI_Method{i}, bestResultsTable.ROI{i}, ...
        bestResultsTable.Best_Correlation(i), bestResultsTable.Best_Concordance(i), bestResultsTable.Best_Time_Interval(i));
end

fclose(fid);
disp('Best LI methods for discordant analyses.');
disp(bestResultsTable);