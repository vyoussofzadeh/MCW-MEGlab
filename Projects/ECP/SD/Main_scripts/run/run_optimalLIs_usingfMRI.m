%%
clc
options = optimset('Display', 'iter', 'TolX', 1e-4);

best_bounds = struct();

ROIs = [1, 2, 6, 11];  % Specifying which ROIs to optimize for
% ROIs = [11];  % Specifying which ROIs to optimize for
windowAdjustments = [0.1:0.01:1.4];  % Example factors to adjust window sizes

% Initialize configuration for each ROI
cfg = struct();
cfg.fmri_LIs_val = fmri_LIs_val;
cfg.wi = wi;
cfg.sub_MF_pt = sub_MF_pt;
cfg.LI_method_label = LI_method_label;
cfg.LI_pt_val_new = LI_pt_val_new;
cfg.fmri_LIs = fmri_LIs;
cfg.IB_megfmri = IB_megfmri;
cfg.rSNR_new = rSNR_new;
cfg.opt_method = opt_method;
cfg.MEG_thre = MEG_thre;
cfg.fMRI_thre = fMRI_thre;
cfg.LI_method_labels = LI_method_labels;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.snr_option = 'thresh'; % 'raw'

% Loop through each ROI
for i = 1:length(ROIs)
    roi = net_sel_mutiple_label{ROIs(i)};  % Extracting the ROI name for indexing
    
    % Initialize best metrics
    bestConcordance = -inf;
    bestBounds = [];
    
    for adjustment = windowAdjustments
        
        cfg.idcx = ROIs(i);
        fun = @(x) - ecpfunc_optimalwindows_dics(setfield(setfield(cfg, 'lowerBound', x(1)), 'upperBound', x(2)));
        
        % Initial guesses and bounds adjusted by current factor
        x0 = [0.4, 1.4] * adjustment;
        lb = [0.2, 0.2] * adjustment;
        ub = [1.8, 1.8] * adjustment;
        
        % Run the optimizer
        [bounds, fval] = fmincon(fun, x0, [], [], [], [], lb, ub, [], options);
        
        cfg_new = cfg; cfg_new.lowerBound = bounds(1); cfg_new.upperBound = bounds(2);
        [MConcordance, bestLIMethod, MCor] = ecpfunc_optimalwindows_dics(cfg_new);
        
        if -fval > bestConcordance
            bestConcordance = -fval;
            bestBounds = bounds;
        end
    end
    
    % Store the best results for current ROI
    best_bounds.(roi) = struct('LowerBound', bestBounds(1), 'UpperBound', bestBounds(2), 'MaxConcordance', bestConcordance, 'BestLIMethod', bestLIMethod, 'BestCorr', MCor);
end

fieldNames = fieldnames(best_bounds);
for i = 1:length(fieldNames)
    roi = fieldNames{i};
    fprintf('ROI: %s\n', roi);
    fprintf('  Lower Bound: %.3f\n', best_bounds.(roi).LowerBound);
    fprintf('  Upper Bound: %.3f\n', best_bounds.(roi).UpperBound);
    fprintf('  Max Concordance: %.3f%%\n', best_bounds.(roi).MaxConcordance);
    fprintf('  Best LI Method: %s\n', best_bounds.(roi).BestLIMethod);
    fprintf('  Best Corr: %.3f\n\n', best_bounds.(roi).BestCorr);
end

fieldNames = fieldnames(best_bounds);  % Get all ROI names
maxConcordances = zeros(length(fieldNames), 1);  % Preallocate array for max concordance values

for i = 1:length(fieldNames)
    roi = fieldNames{i};
    maxConcordances(i) = best_bounds.(roi).MaxConcordance;  % Extract max concordance for each ROI
end

%%
% Plotting
figure;
barHandles = bar(maxConcordances);
title('Maximum Concordance by ROI');
xlabel('ROI');
ylabel('Max Concordance (%)');
set(gca, 'XTickLabel', fieldNames);
xtickangle(45);

% Annotations
xCoords = barHandles(1).XData;
for i = 1:length(maxConcordances)
    text(xCoords(i), maxConcordances(i) + 0.5, sprintf('%.2f%%', maxConcordances(i)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 10, ...
        'Color', 'blue');
    text(xCoords(i), maxConcordances(i) / 2, sprintf('[%.2f, %.2f]', best_bounds.(fieldNames{i}).LowerBound, best_bounds.(fieldNames{i}).UpperBound), ...
        'Rotation', 90, ...
        'HorizontalAlignment', 'right', ...
        'FontSize', 8, ...
        'Color', 'white');
    text(xCoords(i), maxConcordances(i) / 2, sprintf('%s', best_bounds.(fieldNames{i}).BestLIMethod), ...
        'Rotation', 90, ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 8, ...
        'Color', 'green');
end
cfg = []; cfg.outdir = save_dir; filename = 'ConcorMaxOpt'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg); close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');


%%
clc, close all
for methodIdx = 1:length(LI_method_label)
    rSNR_roi = transformPowSubTo3DArrays(rSNR_new.(LI_method_label{methodIdx}));
    
    % Initialize a new matrix to store the absolute differences
    abs_LR = zeros(size(rSNR_roi.left));
    
    % Calculate the absolute difference for each ROI
    for i = 1:size(rSNR_roi.left, 1)  % Iterating over each ROI
        abs_LR(i, :, :) = abs(rSNR_roi.left(i, :, :) - rSNR_roi.right(i, :, :));
        %         abs_LR(i, :, :) = (rSNR_roi.left(i, :, :) );
        
        %     abs_LR(i, :, :) = (rSNR_roi.left(i, :, :) + rSNR_roi.right(i, :, :));
    end
    
    % Assuming abs_LR contains the absolute differences computed as before
    numROIs = size(rSNR_roi.left, 1);
    ROIs = [1,2, 6, 11];
    numSubjects = size(rSNR_roi.left, 2);
    numWindows = size(rSNR_roi.left, 3);
    
    % Matrix to store the 95th percentile for each window in each ROI
    percentile95 = zeros(numROIs, numWindows);
    
    % Calculate the 95th percentile for each window for each ROI
    for i = ROIs
        for j = 1:numWindows
            data_window = squeeze(abs_LR(i, :, j));  % Extract data for all subjects at the current window
            percentile95(i, j) = prctile(data_window, 95);  % Compute the 95th percentile
        end
    end
    
    % Define the threshold as the 95th percentile
    thresholds = 0.2 * max(percentile95, [], 2);  % Calculate max along windows for each ROI
    
    % Finding windows where the mean or median is below the 95th percentile threshold
    best_windows = zeros(numROIs, numWindows);
    for i = ROIs
        for j = 1:numWindows
            % Assuming you are comparing to the mean or median of the absolute differences
            mean_abs_diff = mean(squeeze(abs_LR(i, :, j)));
            if mean_abs_diff > thresholds(i)
                best_windows(i, j) = 1;  % Mark window as "best" if the mean is less than the 95th percentile
            end
        end
    end
    
    % Example: Plotting the 95th percentiles and highlighting the best windows
    best_indices_roi = []; best_bounds  = [];
    
    % close all
    figure;
    j=1;
    for i = ROIs
        subplot(2,2,j), j=j+1;
        roi = net_sel_mutiple_label{i};  % Extracting the ROI name for indexing
        plot(wi(:,1)',percentile95(i, :), 'LineWidth', 2);  % Plot 95th percentile for each ROI
        hold on;
        % Highlight best windows
        best_indices = find(best_windows(i, :) == 1);
        best_indices = [best_indices(1), best_indices(end)];
        %     plot(wi(best_indices,1)', percentile95(i, best_indices), 'ro');
        best_bounds.(roi) = struct('LowerBound', wi(best_indices(1),1), 'UpperBound', wi(best_indices(2),1));
        % Drawing vertical lines for lower and upper bounds
        xline(wi(best_indices(1), 1), 'g', 'LineWidth', 2);  % Green line for lower bound
        xline(wi(best_indices(2), 1), 'b', 'LineWidth', 2);  % Blue line for upper bound
    end
    
    xlabel('Window');
    ylabel('95th Percentile Value');
    title('95th Percentile and Best Windows by ROI');
    legend('95th Percentile', 'Best Windows');
    title(sprintf('Percentile and Bounds for ROI: %s', roi));
    disp(LI_method_label{methodIdx})
    plot_indiv_LI = 0; plot_rSNR = 0; plot_rSNR_LI = 0;
    idcx = [1,2,6,11];
    opt_method = 'rsnr';
    run_optimalLIs_snr_dics_rois
    
    % end
end

%%
plot_indiv_LI = 0; plot_rSNR = 0; plot_rSNR_LI = 0;
idcx = [1,2,6,11];
run_optimalLIs_snr_dics_rois