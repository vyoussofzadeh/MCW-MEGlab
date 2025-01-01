plotCorrMaxByROI(summaryTable, [])

cfg = []; cfg.outdir = save_dir; filename = ' Corr max'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

plotConcMaxByROI(summaryTable, []);
cfg = []; cfg.outdir = save_dir; filename = ' Conc max'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');


%% OLDER VER.
% % Unique ROIs for iteration
% uniqueROIs = unique(summaryTable.ROI);
% 
% % Loop through each ROI for plotting
% figure;
% sgtitle('Corr Max'); % Super title for the figure
% for i = 1:length(uniqueROIs)
%     
%     subplot (4,1,i)
%     roi = uniqueROIs{i};
%     
%     % Extract data for the current ROI
%     roiData = summaryTable(strcmp(summaryTable.ROI, roi), :);
%     
%     % Create figure for current ROI
%     
%     k=1;
%     % Plot Max Values for Correlation
%     hold on;
%     for method = unique(roiData.LI_Method)'
%         methodData = roiData(strcmp(roiData.LI_Method, method) & strcmp(roiData.Metric_Type, 'Correlation'), :);
%         bar(categorical(methodData.LI_Method), methodData.Max_Value, 'BarWidth', 0.2);
%         text(k, methodData.Max_Value + 0.02, sprintf('%.2f', methodData.Max_Value), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
%         k=1+k;
%     end
%     title([roi]);
%     ylabel('Corr.');
%     hold off;
%     set(gca,'color','none');
%     axis tight
%     ylim([0, 1])
%     set(gcf, 'Position', [100, 100, 200, 800]); % Adjust figure size
% end
% 
% cfg = []; cfg.outdir = save_dir; filename = ' Corr max'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
% close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
% 
% 
% figure;
% sgtitle('Conc Max'); % Super title for the figure
% for i = 1:length(uniqueROIs)
%     
%     subplot (4,1, i)
%     % Plot Max Values for Concordance
%     
%     roi = uniqueROIs{i};
%     
%     % Extract data for the current ROI
%     roiData = summaryTable(strcmp(summaryTable.ROI, roi), :);
%     
%     % Create figure for current ROI
%     k=1;
%     hold on;
%     for method = unique(roiData.LI_Method)'
%         methodData = roiData(strcmp(roiData.LI_Method, method) & strcmp(roiData.Metric_Type, 'Concordance'), :);
%         bar(categorical(methodData.LI_Method), methodData.Max_Value, 'BarWidth', 0.2);
%         text(k, methodData.Max_Value + 0.02, sprintf('%.2f', methodData.Max_Value), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
%         k=1+k;
%     end
%     title([roi]);
%     ylabel('Concordance');
%     hold off;
%     set(gca,'color','none');
%     ylim([0, 100])
%     set(gcf, 'Position', [100, 100, 200, 800]); % Adjust figure size    
% end
% 
% cfg = []; cfg.outdir = save_dir; filename = ' Conc max'; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
% close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');
% 
% % disp(summaryTable)