function plotDiscordantTaskPerformance(discordant_subs, meanAccBySubject_Animal, meanAccBySubject_Falsefont, xllabel, ~, ~)
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
discordant_indices_anim = find(ismember(meanAccBySubject_Animal.Subject, discordant_subs_numeric));
discordant_indices_symb = find(ismember(meanAccBySubject_Falsefont.Subject, discordant_subs_numeric));
meanTP_anim = nanmean(meanAccBySubject_Animal.mean_Animal_ACC);
meanTP_symb = nanmean(meanAccBySubject_Falsefont.mean_Falsefont_ACC);

% Plotting Animal Task Performance
figure;
subplot 121
bar(meanAccBySubject_Animal.mean_Animal_ACC(discordant_indices_anim),0.2);
xlabel('Subjects');
set(gca, 'color', 'none');
ylabel('Accuracy (%)');
title({'Animals'});
set(gca, 'XTick', 1:length(discordant_indices_anim), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
hold on; hold off; box off
line(get(gca, 'xlim'), [meanTP_anim meanTP_anim], 'Color', 'red', 'LineStyle', '--');

% Plotting Symbol Task Performance
%     figure;
subplot 122
bar(meanAccBySubject_Falsefont.mean_Falsefont_ACC(discordant_indices_symb),0.2);
xlabel('Subjects');
set(gca, 'color', 'none');
ylabel('Accuracy (%)');
title({'Sybmols'});
set(gca, 'XTick', 1:length(discordant_indices_symb), 'XTickLabel', xllabel, 'XTickLabelRotation', 45);
set(gcf, 'Position', [1000, 100, 600, 300]);
line(get(gca, 'xlim'), [meanTP_symb meanTP_symb], 'Color', 'red', 'LineStyle', '--');
hold on; hold off; box off

end
