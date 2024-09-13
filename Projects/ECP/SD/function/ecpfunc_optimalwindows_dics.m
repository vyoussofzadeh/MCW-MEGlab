function   [MConcordance, BestLIMethod, MCor] = ecpfunc_optimalwindows_dics(cfg_main)

%%
fmri_LIs_val = cfg_main.fmri_LIs_val;
wi = cfg_main.wi;
sub_MF_pt = cfg_main.sub_MF_pt;
LI_method_label = cfg_main.LI_method_label;
LI_pt_val_new = cfg_main.LI_pt_val_new;
fmri_LIs = cfg_main.fmri_LIs;
IB_megfmri = cfg_main.IB_megfmri;
rSNR_new = cfg_main.rSNR_new;
opt_method = cfg_main.opt_method;
lowerBound = cfg_main.lowerBound;
upperBound = cfg_main.upperBound;
MEG_thre = cfg_main.MEG_thre;
fMRI_thre = cfg_main.fMRI_thre;
LI_method_labels = cfg_main.LI_method_labels;
net_sel_mutiple_label = cfg_main.net_sel_mutiple_label;
plot_indiv_LI = cfg_main.plot_indiv_LI;
idcx = cfg_main.idcx;

%% Check if bounds are valid
if cfg_main.lowerBound >= cfg_main.upperBound
    disp('Error: lowerBound must be less than upperBound.');
    MConcordance = 0;  % Assign default or error value to output
    return;  % Exit the function if bounds are not valid
end

%%
% Initialize variables
fMRI_LI = fmri_LIs_val;
timePoints = mean(wi, 2);

sub_IDs = sub_MF_pt;
nsub_IDs = cellfun(@(x) [num2str(find(strcmp(sub_IDs, x))), ':', x], sub_IDs, 'UniformOutput', false);

MEG_LI_ROIs = cell(1, length(idcx));

for j = 1:length(idcx)
    MEG_LI_ROIs{j} = squeeze(LI_pt_val_new.(LI_method_label{1})(idcx(j), :, :)); % Placeholder for methodIdx
end

%%
% Initialize fmri_LIs_ROIs with 11 rows. Assuming the number of columns is known and consistent across the data fields
num_columns = size(fmri_LIs.val.language_Angular, 1);  % Modify as needed based on actual data structure
fmri_LIs_ROIs = NaN(11, num_columns);  % Initialize with NaN for unassigned rows

% Populate specific rows with data
fmri_LIs_ROIs(1,:) = fmri_LIs.val.language_Angular;
fmri_LIs_ROIs(2,:) = fmri_LIs.val.language_Frontal;
fmri_LIs_ROIs(6,:) = fmri_LIs.val.language_Temporal;
fmri_LIs_ROIs(11,:) = fmri_LIs.val.language_Lateral;

% Apply the index mask IB_megfmri to select or reorder rows
% Make sure IB_megfmri contains indices that are valid (1 to 11)
% Assuming IB_megfmri is defined elsewhere and correctly
fmri_LIs_ROIs = fmri_LIs_ROIs(:, IB_megfmri)';  % This line will reorder or filter the rows according to IB_megfmri

%%

summaryTableDynamic = table();

for j = 1:length(idcx)
    
    fMRI_LI = fmri_LIs_ROIs(:, idcx(j));
    
    for methodIdx = 1:length(LI_method_label)
        
        MEG_LI = squeeze(LI_pt_val_new.(LI_method_label{methodIdx})(idcx(j), :, :));
        
        rSNR_roi = transformPowSubTo3DArrays(rSNR_new.(LI_method_label{methodIdx}));
        rSNR_left = squeeze(rSNR_roi.left(idcx(j), :,:));
        rSNR_right = squeeze(rSNR_roi.right(idcx(j), :,:));
        
        if methodIdx == 1, rSNR_left = 20.*rSNR_left; rSNR_right = 20.*rSNR_right; end
        
        switch opt_method
            case 'rsnr'
                [optimalIndices, ~] = findIndividualOptimalTimePoints_interval_rSNR(rSNR_left, rSNR_right, wi, NaN, lowerBound, upperBound);
            case 'LI'
                [optimalIndices, ~] = findIndividualOptimalTimePoints_interval(MEG_LI, fMRI_LI, wi, NaN, lowerBound, upperBound);
        end
        
%         disp(optimalIndices)
%         disp([lowerBound, upperBound])
%         
        optimalInterval = wi(optimalIndices,:);
        meanOptimalTime = mean(optimalInterval);
        
        [concordance, discordantSubs, groupCorrelation, pval] = ...
            calculateConcordanceForTimePoints_interval(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, wi, optimalInterval);
        
        % Store results in the summary table
        newRow = {LI_method_labels{methodIdx}, net_sel_mutiple_label{idcx(j)}, groupCorrelation, concordance, meanOptimalTime, discordantSubs', MEG_LI, fMRI_LI};
        summaryTableDynamic = [summaryTableDynamic; newRow];
        
        if idcx(j) == 11 && plot_indiv_LI == 1 % lateral network (11)
            plotOptimalTimePointsOnMEG2(MEG_LI, fMRI_LI, timePoints,  mean(optimalInterval,2), discordantSubs, MEG_thre, lowerBound, upperBound);
            suptitle(LI_method_label(methodIdx));
            set(gcf, 'Position', [100, 100, 1600, 1300]);
            cfg = []; cfg.outdir = save_dir; filename = ['LI_individuals_',LI_method_label{methodIdx}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
        end
    end
end

% Set column names for the summary table
summaryTableDynamic.Properties.VariableNames = {'LI_Method', 'ROI', 'Correlation', 'Concordance', 'mean_Optimal_Time', 'discord_Subs', 'MEG_LI', 'fMRI_LI'};

summaryTableDynamic_save = summaryTableDynamic;
summaryTableDynamic_save = summaryTableDynamic_save(:,1:end-3);


% Format the numeric data to two decimal places before saving to CSV
summaryTableDynamic_save.Correlation = round(summaryTableDynamic_save.Correlation, 2);
summaryTableDynamic_save.Concordance = round(summaryTableDynamic_save.Concordance, 2);
summaryTableDynamic_save.mean_Optimal_Time = round(summaryTableDynamic_save.mean_Optimal_Time, 2);


%% MAX conc
% Determine the maximum concordance and the method that produced it
[MConcordance, maxIdx] = max(summaryTableDynamic.Concordance);
BestLIMethod = summaryTableDynamic.LI_Method{maxIdx};
% MCor = summaryTableDynamic.Correlation(maxIdx);
[MCor, ~] = max(summaryTableDynamic.Correlation);
% MConcordance = max(summaryTableDynamic_save.Correlation);


%% MAX Corr
% [MConcordance, maxIdx] = max(summaryTableDynamic.Correlation);
% BestLIMethod = summaryTableDynamic.LI_Method{maxIdx};
% % MCor = summaryTableDynamic.Correlation(maxIdx);
% [MCor, ~] = max(summaryTableDynamic.Concordance);

end