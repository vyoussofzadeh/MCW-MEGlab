% ------------------------------------------------------------------------
metricNames       = {'Correlation','Concordance'};
roi_labels        = {'Ang','Front','Temp','Lat'};
network_sel       = 1:numel(roi_labels);                 % ROIs to plot
LI_method_labels  = {'SourceMag','Count','Bootstrp'};
midpoints         = mean(wi,2);                          % time axis

customColors = [  0 114 189 ;
    217  83  25 ;
    237 177  32 ;
    126  47 142 ]/256;                       % 1 colour / ROI

close all
for metricIdx = 1:numel(metricNames)
    figure
    hold on
    legendEntries = cell(numel(network_sel),1);
    
    for net_idx = 1:numel(network_sel)
        roiID = network_sel(net_idx);
        
        % -------- gather 3×time matrix for this ROI ---------------------
        LI_values = nan(numel(LI_method_labels), numel(midpoints));
        for m = 1:numel(LI_method_labels)
            M = resultsTable.Metrics(m).(metricNames{metricIdx});  % [ROI × time]
            LI_values(m,:) = M(roiID,:);
        end
        
        % -------- shaded mean ± SEM -------------------------------------
        AlphaLine(midpoints, LI_values, customColors(net_idx,:), 'LineWidth',1.6);
        
        % -------- find & annotate peak of mean curve --------------------
        mu = mean(LI_values,1,'omitnan');
        [~,pk] = max(mu);
        pkTime = midpoints(pk);
        pkVal  = mu(pk);
        
        % build the annotation string once per metric ---------------------------
        switch metricNames{metricIdx}
            case 'Correlation'          % two decimals, keep sign
                lbl = sprintf('%.2f @ %.1fs', pkVal, pkTime);
            case 'Concordance'          % whole-number percent
                lbl = sprintf('%.0f%% @ %.1fs', pkVal, pkTime);
        end
        
        text(pkTime, pkVal, lbl, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom');
        
        
        % label
        %         text(pkTime, pkVal, sprintf('%.1fs', pkTime), ...
        %             'Horiz','center','Vert','bottom');
        %
        % *** vertical reference line (new) ***
%         line([pkTime pkTime], ylim, ...
%             'Color', customColors(net_idx,:), ...
%             'LineWidth',1.3, 'LineStyle','--');
        
        hRef = line([pkTime pkTime], ylim, ...
            'Color', customColors(net_idx,:), ...
            'LineWidth',1.3,'LineStyle','--', ...
            'HandleVisibility','off');    % <-- legend ignores this
        
        legendEntries{net_idx} = roi_labels{roiID};
    end
    
    % -------- axis cosmetics -------------------------------------------
    xlabel('time (s)')
    switch metricNames{metricIdx}
        case 'Correlation',  ylabel('r'); ylim([-0.3 1])
        case 'Concordance',  ylabel('%'); ylim([0 90])
    end
    title(metricNames{metricIdx})
%     box off
    set(gca, 'color', 'none');
    
    
    legend(legendEntries,'Location','southoutside', ...
        'Orientation','horizontal','NumColumns',numel(network_sel))
    box off;
    set(gcf, 'Position', [800, 400, 500, 800]);
    hold off; % Release the plot hold
    axis tight
    
    xlabel('Time (s)');
    set(gca, 'color', 'none');
    
    % --------------- export this figure --------------
     
    cfg = [];
    cfg.outdir   = save_dir;
    cfg.filename = sprintf('MEG_fMRI_LI_%s_shadedPerROI_vertical', ...
        lower(metricNames{metricIdx}));  % unique name
    cfg.type     = 'svg';
    do_export_fig(cfg);
    close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
    
end
