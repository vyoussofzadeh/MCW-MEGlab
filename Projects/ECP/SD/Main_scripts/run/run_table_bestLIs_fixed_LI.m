% Within your existing for loop over ROIs and methods:

summaryTableFixed = table();


for j = 1:length(network_sel)
    
    fMRI_LI = fmri_LIs_ROIs(:, j);

    for methodIdx = 1:length(LI_method_label)
        
        MEG_LI = squeeze(LI_pt_val_new.(LI_method_label{methodIdx})(network_sel(j), :, :));
        
        rSNR_roi = transformPowSubTo3DArrays(rSNR_new.(LI_method_label{methodIdx}));
        rSNR_left = squeeze(rSNR_roi.left(network_sel(j), :,:));
        rSNR_right = squeeze(rSNR_roi.right(network_sel(j), :,:));
        
        if methodIdx == 1
            rSNR_left = 20 .* rSNR_left;
            rSNR_right = 20 .* rSNR_right;
        end

        roiName = net_sel_mutiple_label{network_sel(j)};

        switch opt_method
            case 'fixed'
                % 1) Extract the columns of MEG_LI for each row in summaryTable 
                %    that matches these "fixed" intervals
                MEG_LI_cell = getMEGLIForIntervals(MEG_LI, summaryTable, wi);

                % 2) Optionally, loop over each row in summaryTable or 
                %    just pick the row(s) relevant to this ROI & method

                % Find rows of summaryTable that match the current method & ROI, e.g.:
                rowsMask = strcmp(summaryTable.LI_Method, LI_method_label{methodIdx}) ...
                         & strcmp(summaryTable.ROI, roiName) & strcmp(summaryTable.Metric_Type, 'Concordance');

                theseRows = find(rowsMask);
                % For each row in summaryTable that matches ROI & method:
                for rIdx = theseRows'
                    
                    % The sub-matrix for that row (all subjects x columns for that interval)
                    MEG_LI_sub = MEG_LI_cell{rIdx};
                    % If empty, means no exact match found in getMEGLIForIntervals
                    if isempty(MEG_LI_sub)
                        disp('No matching interval in wi for:');
                        disp(summaryTable.Time_Interval(rIdx,:));
                        continue;
                    end

                    % Suppose you want the time interval from summaryTable:
                    fixedInterval = summaryTable.Time_Interval(rIdx,:);
 
                    % Now you can run any calculations using that interval
                    % For instance, if you want to treat all columns in MEG_LI_sub as "valid",
                    % define an interval in 'wi' or skip the function that needs intervals.
                    
                    % Example: we can compute a *single* average across those columns 
                    % for each subject:
                    fixedMEG_LI = mean(MEG_LI_sub, 2, 'omitnan');  % (#Subjects x 1)

                    % Then compute some measure of correlation or 'concordance' with fMRI_LI
                    % For your existing code, you might do something like:
                    [concordance, discordantSubs, groupCorrelation, ...
                     optimalMEG_LI, pval, kappa] = calculateConcordanceForTimePoints_interval( ...
                           MEG_LI, MEG_thre, fMRI_LI, fMRI_thre, wi, fixedInterval);
                       
                       
                    % Store or do something with the results...
                    newRow = {
                        LI_method_label{methodIdx}, ...
                        roiName, ...
                        groupCorrelation, ...
                        pval, ...
                        concordance, ...
                        kappa, ...
                        fixedInterval, ...
                        discordantSubs', ...
                        MEG_LI, ...
                        optimalMEG_LI, ...
                        fMRI_LI, ...
                        rSNR_left, ...
                        rSNR_right };

                    summaryTableFixed = [summaryTableFixed; newRow];
                end
        end
    end
end

summaryTableFixed.Properties.VariableNames = {'LI_Method', 'ROI', 'Correlation', 'Corr_P_value', 'Concordance', 'kappa', 'mean_Optimal_Time', 'discord_Subs', 'MEG_LI', 'optimalMEG_LI', 'fMRI_LI', 'rSNR_left', 'rSNR_right'};

disp('Summary table for fixed intervals.');
disp(summaryTableFixed)


