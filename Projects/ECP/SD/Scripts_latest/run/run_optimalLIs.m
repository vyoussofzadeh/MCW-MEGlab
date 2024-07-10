close all,

% Initialize variables
fMRI_LI = fmri_LIs_val;
timePoints = mean(wi, 2);
% lowerBound = 0.45; upperBound = 1.1;
% lowerBound = 0.4; upperBound = 0.9;
% lowerBound = 0.4; upperBound = 0.9;
% lowerBound = 0.2; upperBound = 1.1;
% lowerBound = 0.35; upperBound = 1.0;
% lowerBound = 0.4; upperBound = 1.0;
% lowerBound = 0.3; upperBound = 1.4;

sub_IDs = sub_MF_pt;
nsub_IDs = cellfun(@(x) [num2str(find(strcmp(sub_IDs, x))), ':', x], sub_IDs, 'UniformOutput', false);

idcx = [1, 2, 6, 11];
MEG_LI_ROIs = cell(1, length(idcx));

for j = 1:length(idcx)
    MEG_LI_ROIs{j} = squeeze(LI_pt_val_new.(LI_method_label{1})(idcx(j), :, :)); % Placeholder for methodIdx
end

fmri_LIs_ROIs = [fmri_LIs.val.language_Angular, fmri_LIs.val.language_Frontal, ...
    fmri_LIs.val.language_Temporal, fmri_LIs.val.language_Lateral];
fmri_LIs_ROIs = fmri_LIs_ROIs(IB_megfmri,:);

summaryTableDynamic = table();

for j = 1:length(idcx)
    
    fMRI_LI = fmri_LIs_ROIs(:, j);
    
    for methodIdx = 1:length(LI_method_label)
        
        MEG_LI = squeeze(LI_pt_val_new.(LI_method_label{methodIdx})(idcx(j), :, :));
        
        % Compute group-level correlation
        [groupCorrelation, optimalInterval, optimalInterval_all, pval]= ...
            computeGroupLevelMEGfMRICorrelation_timepoints_interval(MEG_LI, fMRI_LI, wi, lowerBound, upperBound);
        disp(['pval:', num2str(pval)])
        
        meanOptimalTime = mean(optimalInterval);
        
        [concordance, discordantSubs, ~] = ...
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

% Save summary table
writetable(summaryTableDynamic, 'LI_Metrics_Summary_Dynamic.csv');

% Save dynamic summary table as text file
fid = fopen('LI_Metrics_Summary_Dynamic.txt', 'wt');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', summaryTableDynamic.Properties.VariableNames{:});

for i = 1:height(summaryTableDynamic)
    fprintf(fid, '%s\t%s\t%f\t%f\t%f\n', summaryTableDynamic.LI_Method{i}, summaryTableDynamic.ROI{i}, ...
        summaryTableDynamic.Correlation(i), summaryTableDynamic.Concordance(i), summaryTableDynamic.mean_Optimal_Time(i));
end

fclose(fid);
disp('Summary table for dynamic intervals.');
disp(summaryTableDynamic)