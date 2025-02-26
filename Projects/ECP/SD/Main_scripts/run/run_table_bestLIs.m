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
            Corr_P_value = methodData.Corr_P_value;
            bestMethod = methodData.LI_Method;
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
    newRow = {bestMethod, roi, bestCorrelation, Corr_P_value, bestConcordance, bestTimeInterval, bestDiscordSubs, bestDiscordSubs_megli, bestDiscordSubs_optmegli, bestDiscordSubs_fmrili, bestDiscordSubs_rSNR_L, bestDiscordSubs_rSNR_R};
    bestResultsTable = [bestResultsTable; newRow];
end

% Set column names for the best results table
% bestResultsTable.Properties.VariableNames = {'Best_LI_Method', 'ROI', 'Best_Correlation', 'Best_Concordance', 'Best_Time_Interval', 'Best_Discord_Subs', 'MEG_LI', 'optMEG_LI', 'fMRI_LI', 'rSNR_L', 'rSNR_R'};
bestResultsTable.Properties.VariableNames = {'Best_LI_Method', 'ROI', 'Best_Correlation', 'Corr_pval', 'Best_Concordance', 'Best_Time_Interval', 'Best_Discord_Subs', 'MEG_LI', 'optMEG_LI', 'fMRI_LI', 'rSNR_L', 'rSNR_R'};

bestResultsTable_save = bestResultsTable;
bestResultsTable_save = bestResultsTable_save(:,1:end-3);

% Format the numeric data to two decimal places before saving to CSV
bestResultsTable_save.Best_Correlation = round(bestResultsTable_save.Best_Correlation, 2);
bestResultsTable_save.Best_Concordance = round(bestResultsTable_save.Best_Concordance, 2);
bestResultsTable_save.Best_Time_Interval = round(bestResultsTable_save.Best_Time_Interval, 2);

% Save best results table
writetable(bestResultsTable_save, fullfile(save_dir,'Best_LI_Methods_Summary.csv'));


% Save best results table as text file
fid = fopen(fullfile(save_dir,'Best_LI_Methods_Summary.txt'), 'wt');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', bestResultsTable.Properties.VariableNames{:});

for i = 1:height(bestResultsTable)
    fprintf(fid, '%s\t%s\t%.2f\t%.2f\t%.2f\n', bestResultsTable.Best_LI_Method{i}, bestResultsTable.ROI{i}, ...
        bestResultsTable.Best_Correlation(i), bestResultsTable.Best_Concordance(i), bestResultsTable.Best_Time_Interval(i));
end

fclose(fid);

%% grossDiscordSubs (left and right discordance)
clc
% Create a new column in bestResultsTable to hold subject indices
bestResultsTable.Gross_Discord_Subs = cell(height(bestResultsTable),1);
bestResultsTable.MEG_LI_tern  = cell(height(bestResultsTable),1);
bestResultsTable.fMRI_LI_tern = cell(height(bestResultsTable),1);

for iRow = 1:height(bestResultsTable)
    % Extract the MEG & fMRI LI arrays (Nx1 each)
    optimalMEG_LI = bestResultsTable.optMEG_LI{iRow};
    fmri_LI       = bestResultsTable.fMRI_LI{iRow};
    
    % 1) Ternary classification for MEG
    cfg = [];
    cfg.thre = MEG_thre;          % e.g. 10
    cfg.LI   = optimalMEG_LI;     % Nx1 numeric
    MEG_LI_tern = do_ternary_classification2(cfg);
    % e.g. MEG_LI_tern is Nx1 numeric with 0=Right, 1=Left, 2=Mid/Unknown
    
    % 2) Ternary classification for fMRI
    cfg = [];
    cfg.thre = fMRI_thre;         % e.g. 10
    cfg.LI   = fmri_LI;
    fMRI_LI_tern = do_ternary_classification2(cfg);
    % Nx1 numeric (0=Right, 1=Left, 2=Mid)
    
    % --- Store these ternary vectors in the new columns
    bestResultsTable.MEG_LI_tern{iRow}  = MEG_LI_tern;
    bestResultsTable.fMRI_LI_tern{iRow} = fMRI_LI_tern;
    
    grossDiscordSubs = find( (MEG_LI_tern==1 & fMRI_LI_tern==-1) | ...
        (MEG_LI_tern==-1 & fMRI_LI_tern==1) );
    
    % Store these subject indices in the new column
    bestResultsTable.Gross_Discord_Subs{iRow} = grossDiscordSubs;
    
    % (Optional) Print them
    %     fprintf('ROI #%d = %s: #GrossDiscord = %d\n', iRow, bestResultsTable.ROI{iRow}, length(grossDiscordSubs));
    %     disp(grossDiscordSubs);
    
end

%%
% % 1) Initialize the columns as cell arrays
% %    (One cell per row, each containing the Nx1 ternary vector)
% bestResultsTable.MEG_LI_tern  = cell(height(bestResultsTable),1);
% bestResultsTable.fMRI_LI_tern = cell(height(bestResultsTable),1);
% 
% % 2) For loop across each ROI row in bestResultsTable
% for iRow = 1:height(bestResultsTable)
%     % Extract the MEG & fMRI LI arrays (Nx1 each)
%     optimalMEG_LI = bestResultsTable.optMEG_LI{iRow};  % numeric Nx1
%     fmri_LI       = bestResultsTable.fMRI_LI{iRow};   % numeric Nx1
%     
%     % --- MEG ternary classification ---
%     cfg = [];
%     cfg.thre = MEG_thre;  % e.g. 10
%     cfg.LI   = optimalMEG_LI;
%     MEG_LI_tern = do_ternary_classification2(cfg); 
%     % MEG_LI_tern is Nx1 numeric, e.g. -1=Right, +1=Left, 0=Mid
%     
%     % --- fMRI ternary classification ---
%     cfg = [];
%     cfg.thre = fMRI_thre; % e.g. 10
%     cfg.LI   = fmri_LI;
%     fMRI_LI_tern = do_ternary_classification2(cfg);
% 
%     % --- Store these ternary vectors in the new columns
%     bestResultsTable.MEG_LI_tern{iRow}  = MEG_LI_tern;
%     bestResultsTable.fMRI_LI_tern{iRow} = fMRI_LI_tern;
%     
%     % --- (Optional) compute 'Gross_Discord_Subs' based on Left vs. Right mismatch
%     grossDiscordSubs = find( ...
%         (MEG_LI_tern == 1 & fMRI_LI_tern == -1) | ...
%         (MEG_LI_tern == -1 & fMRI_LI_tern == 1) ...
%     );
%     bestResultsTable.Gross_Discord_Subs{iRow} = grossDiscordSubs;
% end

%% Gross discordance (blandAltman)
% bestResultsTable.Gross_Discord_Subs = cell(height(bestResultsTable),1);
% 
% for iRow = 1:height(bestResultsTable)
%     % Extract the MEG & fMRI LI arrays
%     optimalMEG_LI = bestResultsTable.optMEG_LI{iRow};  % e.g. Nx1
%     fmri_LI       = bestResultsTable.fMRI_LI{iRow};   % e.g. Nx1
%     
%     
%     cfg = []; cfg.thre = MEG_thre;
%     cfg.LI = optimalMEG_LI; MEG_LI_tern = do_ternary_classification2(cfg);
%     
%     cfg = []; cfg.thre = fMRI_thre;
%     cfg.LI = fmri_LI; fMRI_LI_tern = do_ternary_classification2(cfg);
% 
%     
%     % If they are bigger (NxT), pick a column or average, etc.
%     % But assuming Nx1:
%     [concordanceBA, outlierIdx] = ecpfunc_blandAltmanConcordance(optimalMEG_LI, fmri_LI);
%     
%     % Store the outlier indices in the new column:
%     bestResultsTable.Gross_Discord_Subs{iRow} = outlierIdx;
%     
%     % (Optional) display or store the Concordance, too, if you want
%     % bestResultsTable.BA_Concordance(iRow) = concordanceBA;
%     
%     fprintf('ROI=%s, Method=%s | Concordance=%.2f, #GrossOutliers=%d\n', ...
%         bestResultsTable.ROI{iRow}, bestResultsTable.Best_LI_Method{iRow}, ...
%         100*concordanceBA, length(outlierIdx));
% end

%%
disp('Best LI methods for discordant analyses.');
disp(bestResultsTable);