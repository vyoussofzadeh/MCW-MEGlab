function plotDiscordantTaskPerformance2(discordant_subs, meanAccBySubject_Animal, meanAccBySubject_Falsefont, xllabel, save_dir, method)
% plotDiscordantTaskPerformance: Plot task performance for discordant subjects
%
% Inputs:
%   - discordant_subs: Cell array of discordant subject IDs as strings
%   - meanAccBySubject_Animal: Table containing accuracy data for the animal task
%   - meanAccBySubject_Falsefont: Table containing accuracy data for the symbol task
%   - xllabel: Cell array of labels for the x-axis
%   - save_dir: Directory to save the plot
%   - LI_method_label: Cell array of LI method labels
%   - LI_method: Index of the current LI method

% Convert discordant subjects to numeric
discordant_subs_numeric = cellfun(@(x) str2double(x(3:end)), discordant_subs);
discordant_indices_anim = ismember(meanAccBySubject_Animal.Subject, discordant_subs_numeric);
discordant_indices_symb = ismember(meanAccBySubject_Falsefont.Subject, discordant_subs_numeric);
meanTP_anim = nanmean(meanAccBySubject_Animal.mean_Animal_ACC);
meanTP_symb = nanmean(meanAccBySubject_Falsefont.mean_Falsefont_ACC);

% % Plotting Animal Task Performance
% figure;
% subplot 121
% bar(meanAccBySubject_Animal.mean_Animal_ACC(discordant_indices_anim));
% xlabel('Subjects');
% set(gca, 'color', 'none');
% ylabel('Accuracy (%)');
% title({'Animals'});
% set(gca, 'XTick', 1:length(discordant_indices_anim), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
% %     set(gcf, 'Position', [1000, 100, 300, 300]);
% hold on;
% line(get(gca, 'xlim'), [meanTP_anim meanTP_anim], 'Color', 'red', 'LineStyle', '--');
% hold off;
% cfg = [];
% cfg.outdir = save_dir;
% filename = ['TaskPerformace_anim_Dicordant_LIs_', method];
% cfg.filename = filename;
% cfg.type = 'fig';
% do_export_fig(cfg);
% 
% % Plotting Symbol Task Performance
% %     figure;
% subplot 122
% bar(meanAccBySubject_Falsefont.mean_Falsefont_ACC(discordant_indices_symb));
% xlabel('Subjects');
% set(gca, 'color', 'none');
% ylabel('Accuracy (%)');
% title({'Sybmols'});
% set(gca, 'XTick', 1:length(discordant_indices_symb), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
% set(gcf, 'Position', [1000, 100, 600, 300]);
% hold on;
% line(get(gca, 'xlim'), [meanTP_symb meanTP_symb], 'Color', 'red', 'LineStyle', '--');
% hold off;


%     cfg = [];
%     cfg.outdir = save_dir;
%     filename = ['TaskPerformace_symb_Dicordant_LIs_', method];
%     cfg.filename = filename;
%     cfg.type = 'fig';
%     do_export_fig(cfg);

%%

figure;
subplot 121
daboxplot(meanAccBySubject_Animal.mean_Animal_ACC, 'groups', ones(1, numel(meanAccBySubject_Animal.mean_Animal_ACC)), 'outsymbol', 'kx', 'xtlabels', 'test', 'fill', 0);
xlabel('All Subjects');
ylabel('Accuracy (%)');
set(gca, 'FontSize', 10);
discordant_ACC_anim = meanAccBySubject_Animal.mean_Animal_ACC(discordant_indices_anim);
hold on;
hScatter = scatter(ones(size(discordant_ACC_anim)), discordant_ACC_anim, 'r', 'filled');
hold off;
title({'Animals'});
l = legend(hScatter, 'Discordant', 'Location', 'southout');
% legendPos = l.Position;
% legendPos(1) = legendPos(1) + 0.02;
% l.Position = legendPos;
set(gca, 'XTick', []);
box off;
% set(gcf, 'Position', [1000   100   300   300]);
set(gca,'color','none');

% cfg = []; cfg.outdir = save_dir; filename = ['2_TaskPerformace_anim_Dicordant_LIs_', LI_method_label{LI_method}]; cfg.filename = filename; cfg.type = 'fig'; do_export_fig(cfg)

% figure;
subplot 122
daboxplot(meanAccBySubject_Falsefont.mean_Falsefont_ACC, 'groups', ones(1, numel(meanAccBySubject_Falsefont.mean_Falsefont_ACC)), 'outsymbol', 'kx', 'xtlabels', 'test', 'fill', 0);
xlabel('All Subjects');
ylabel('Accuracy (%)');
set(gca, 'FontSize', 10);
discordant_ACC_symb = meanAccBySubject_Falsefont.mean_Falsefont_ACC(discordant_indices_symb);
hold on;
hScatter = scatter(ones(size(discordant_ACC_symb)), discordant_ACC_symb, 'r', 'filled');
hold off;
title({'Sybmols'});
l = legend(hScatter, 'Discordant', 'Location', 'southout');
% legendPos = l.Position;
% legendPos(1) = legendPos(1) + 0.02;
% l.Position = legendPos;
set(gca, 'XTick', []);
box off;
set(gcf, 'Position', [1000, 100, 600, 300]);
set(gca,'color','none');

end
