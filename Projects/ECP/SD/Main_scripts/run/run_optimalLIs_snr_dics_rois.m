
% Initialize variables
fMRI_LI = fmri_LIs_val;
timePoints = wi(:,1);

sub_IDs = sub_MF_pt;
nsub_IDs = cellfun(@(x) [num2str(find(strcmp(sub_IDs, x))), ':', x], sub_IDs, 'UniformOutput', false);

MEG_LI_ROIs = cell(1, length(network_sel));

for j = 1:length(network_sel)
    MEG_LI_ROIs{j} = squeeze(LI_pt_val_new.(LI_method_label{1})(network_sel(j), :, :)); % Placeholder for methodIdx
end

fmri_LIs_ROIs = [fmri_LIs.val.language_Angular, fmri_LIs.val.language_Frontal, ...
    fmri_LIs.val.language_Temporal, fmri_LIs.val.language_Lateral];
fmri_LIs_ROIs = fmri_LIs_ROIs(IB_megfmri,:);

summaryTableDynamic = table();

for j = 1:length(network_sel)
    
    fMRI_LI = fmri_LIs_ROIs(:, j);
    
    for methodIdx = 1:length(LI_method_label)
        
        MEG_LI = squeeze(LI_pt_val_new.(LI_method_label{methodIdx})(network_sel(j), :, :));
        
        rSNR_roi = transformPowSubTo3DArrays(rSNR_new.(LI_method_label{methodIdx}));
        rSNR_left = squeeze(rSNR_roi.left(network_sel(j), :,:));
        rSNR_right = squeeze(rSNR_roi.right(network_sel(j), :,:));
        
        if methodIdx == 1, rSNR_left = 20.*rSNR_left; rSNR_right = 20.*rSNR_right; end
        
        plot_flag = 1;
        
        roiName = net_sel_mutiple_label{network_sel(j)};
        
        switch opt_method
            case 'rsnr'
                lowerBound = best_bounds.(roiName).LowerBound; upperBound = best_bounds.(roiName).UpperBound;
                [optimalIndices, maxDiff_opt] = findIndividualOptimalTimePoints_interval_rSNR(rSNR_left, rSNR_right, wi, NaN, lowerBound, upperBound);
            case 'LI'
                lowerBound = best_bounds.(roiName).LowerBound; upperBound = best_bounds.(roiName).UpperBound;
                [optimalIndices, maxDiff_opt] = findIndividualOptimalTimePoints_interval(MEG_LI, fMRI_LI, wi, NaN, lowerBound, upperBound);
            case 'AUC'
%                 [optimalIndices, maxAUC_opt] = findIndividualOptimalTimePoints_maxAUC(rSNR_left, rSNR_right, wi, NaN, lowerBound, upperBound);
                [optimalIndices, maxDiff_opt, bounds] = findIndividualOptimalTimePoints_interval_rSNR_slidingWindow(rSNR_left, rSNR_right, wi, NaN, minlowerband, maxUpperband);
            case 'DomH'
                [optimalIndices, maxSNR_opt] = findIndividualOptimalTimePoints_dominantHemisphere(rSNR_left, rSNR_right, wi, NaN, lowerBound, upperBound);
            case 'rsnrmax'
                [optimalIndices, maxAUC_opt] = findIndividualOptimalTimePoints_interval_rSNR_MaxInd(rSNR_left, rSNR_right, wi, NaN, lowerBound, upperBound);
            case 'wrsnr'
                [optimalIndices, maxDiff_opt] = findIndividualOptimalTimePoints_interval_WrSNR(rSNR_left, rSNR_right, wi, NaN, lowerBound, upperBound);
            case {'rsnr_optbound', 'rsnr_optbound_mean'}
                [optimalIndices, maxDiff_opt, bounds] = findIndividualOptimalTimePoints_interval_rSNR_optbound(rSNR_left, rSNR_right, wi, NaN, minlowerband, maxUpperband);
        end
        
        optimalInterval = wi(optimalIndices,:);
        meanOptimalTime = mean(optimalInterval);
        
        [concordance, discordantSubs, groupCorrelation, optimalMEG_LI, pval, kappa] = ...
            calculateConcordanceForTimePoints_interval(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, wi, optimalInterval);
        
        switch opt_method
            case 'rsnr_optbound_mean'
                optimalInterval = wi(bounds,:);
                [concordance, discordantSubs, groupCorrelation, optimalMEG_LI, pval] = ...
                    calculateConcordanceForTimePoints_interval(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, wi, optimalInterval);
        end
        
        % Store results in the summary table
        newRow = {LI_method_labels{methodIdx}, net_sel_mutiple_label{network_sel(j)}, groupCorrelation, pval, concordance, kappa, meanOptimalTime, discordantSubs', MEG_LI, optimalMEG_LI, fMRI_LI, rSNR_left, rSNR_right};      
        summaryTableDynamic = [summaryTableDynamic; newRow];
        
        if network_sel(j) == 11 && plot_indiv_LI == 1 % lateral network (11)
            plotOptimalTimePointsOnMEG2(MEG_LI, fMRI_LI, timePoints,  mean(optimalInterval,2), discordantSubs, MEG_thre, lowerBound, upperBound);
            suptitle(LI_method_label(methodIdx));
            set(gcf, 'Position', [100, 100, 1600, 1300]);
            cfg = []; cfg.outdir = save_dir; filename = ['LI_individuals_',LI_method_label{methodIdx}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
        end
    end
end

% Set column names for the summary table
summaryTableDynamic.Properties.VariableNames = {'LI_Method', 'ROI', 'Correlation', 'Corr_P_value', 'Concordance', 'kappa', 'mean_Optimal_Time', 'discord_Subs', 'MEG_LI', 'optimalMEG_LI', 'fMRI_LI', 'rSNR_left', 'rSNR_right'};
summaryTableDynamic.Corr_P_value = round(summaryTableDynamic.Corr_P_value, 6);
summaryTableDynamic_save = summaryTableDynamic;
summaryTableDynamic_save = summaryTableDynamic_save(:,1:end-5);


% Format the numeric data to two decimal places before saving to CSV
summaryTableDynamic_save.Correlation = round(summaryTableDynamic_save.Correlation, 2);
summaryTableDynamic_save.Concordance = round(summaryTableDynamic_save.Concordance, 2);
summaryTableDynamic_save.mean_Optimal_Time = round(summaryTableDynamic_save.mean_Optimal_Time, 2);

% Save summary table
writetable(summaryTableDynamic_save, fullfile(save_dir,'LI_Metrics_Summary_Dynamic.csv'));

% % Save dynamic summary table as text file
fid = fopen(fullfile(save_dir,'LI_Metrics_Summary_Dynamic.txt'), 'wt');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', summaryTableDynamic_save.Properties.VariableNames{:});

for i = 1:height(summaryTableDynamic_save)
    fprintf(fid, '%s\t%s\t%f\t%f\t%f\n', summaryTableDynamic_save.LI_Method{i}, summaryTableDynamic_save.ROI{i}, ...
        summaryTableDynamic_save.Correlation(i), summaryTableDynamic_save.Concordance(i), summaryTableDynamic_save.mean_Optimal_Time(i));
end

fclose(fid);
disp('Summary table for dynamic intervals.');
disp(summaryTableDynamic_save)

%% Plot rSNR
if plot_rSNR == 1
    for j = 4:4 %length(network_sel)
        fMRI_LI = fmri_LIs_ROIs(:, j);
        for methodIdx = 1:length(LI_method_label)
            MEG_LI = squeeze(LI_pt_val_new.(LI_method_label{methodIdx})(network_sel(j), :, :));
            rSNR_roi = transformPowSubTo3DArrays(rSNR_new.(LI_method_label{methodIdx}));
            rSNR_left = squeeze(rSNR_roi.left(network_sel(j), :,:));
            rSNR_right = squeeze(rSNR_roi.right(network_sel(j), :,:));
            if methodIdx == 1, rSNR_left = 20.*rSNR_left; rSNR_right = 20.*rSNR_right; end
            
            roiName = net_sel_mutiple_label{network_sel(j)};
            
            
            switch opt_method
                case {'rsnr', 'LI'}
                    lowerBound = best_bounds.(roiName).LowerBound; upperBound = best_bounds.(roiName).UpperBound;
                    [optimalIndices, maxDiff_opt] = findIndividualOptimalTimePoints_interval_rSNR3(rSNR_left, rSNR_right, fMRI_LI, wi, 1:length(rSNR_right), lowerBound, upperBound);
                case 'rsnr_optbound'
                    [optimalIndices, maxDiff_opt] = findIndividualOptimalTimePoints_interval_rSNR4(rSNR_left, rSNR_right, fMRI_LI, wi, 1:length(rSNR_right), bounds);
            end
            
            sgtitle(LI_method_label(methodIdx))
            cfg = []; cfg.outdir = save_dir; filename = ['rSNR_individuals_',LI_method_label{methodIdx}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
        end
    end
end


%% Plot rSNR+LI
if plot_rSNR_LI == 1
    for j = 4:4 %length(network_sel)
        fMRI_LI = fmri_LIs_ROIs(:, j);
        for methodIdx = 1:length(LI_method_label)
            
            MEG_LI = squeeze(LI_pt_val_new.(LI_method_label{methodIdx})(network_sel(j), :, :));
            rSNR_roi = transformPowSubTo3DArrays(rSNR_new.(LI_method_label{methodIdx}));
            rSNR_left = squeeze(rSNR_roi.left(network_sel(j), :,:));
            rSNR_right = squeeze(rSNR_roi.right(network_sel(j), :,:));
            if methodIdx == 1, rSNR_left = 20.*rSNR_left; rSNR_right = 20.*rSNR_right; end
            
            roiName = net_sel_mutiple_label{network_sel(j)};
            
            switch opt_method
                case {'rsnr','LI'}
                    lowerBound = best_bounds.(roiName).LowerBound; upperBound = best_bounds.(roiName).UpperBound;
                    [optimalIndices, maxDiff_opt] = findIndividualOptimalTimePoints_interval_rSNR(rSNR_left, rSNR_right, wi, NaN, lowerBound, upperBound);
                case 'rsnr_optbound'
                    [optimalIndices, maxDiff_opt, bounds] = findIndividualOptimalTimePoints_interval_rSNR_optbound(rSNR_left, rSNR_right, wi, NaN, minlowerband, maxUpperband);
            end
            
            optimalInterval = wi(optimalIndices,:);
            
            [concordance, discordantSubs, groupCorrelation, pval] = calculateConcordanceForTimePoints_interval(MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, wi, optimalInterval);
            
            rSNR_MEG = []; rSNR_MEG.rSNR_left = rSNR_left; rSNR_MEG.rSNR_right = rSNR_right;
            
            plotOptimalTimePointsOnMEG4(rSNR_MEG, MEG_LI, fMRI_LI, wi, optimalIndices, discordantSubs, MEG_thre, bounds); sgtitle(LI_method_label(methodIdx))
%             plotOptimalTimePointsOnMEG4(rSNR_MEG, MEG_LI(72,:), fMRI_LI(72), wi, optimalIndices, discordantSubs, MEG_thre, bounds(72,:)); 
%             sgtitle(LI_method_label(methodIdx))

%             cfg = []; cfg.outdir = save_dir; filename = ['LI_rSNR_individuals_',LI_method_label{methodIdx}]; cfg.filename = filename; cfg.type = 'png'; do_export_fig(cfg); combined_path = fullfile(save_dir,[cfg.filename, '.png']); 
            cfg = []; cfg.outdir = save_dir; filename = ['LI_rSNR_individuals_',LI_method_label{methodIdx}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

            
            plotOptimalTimePointsOnMEG4_selective(rSNR_MEG, MEG_LI, fMRI_LI, wi, optimalIndices, discordantSubs, MEG_thre, bounds);
            sgtitle(LI_method_label(methodIdx))
            cfg = []; cfg.outdir = save_dir; filename = ['LI_rSNR_discordance_',LI_method_label{methodIdx}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
            
        end
    end
end