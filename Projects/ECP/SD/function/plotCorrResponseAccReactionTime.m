function plotCorrResponseAccReactionTime(accuracyResults_updt, T_patn_MEGfMRI, discordSubs)
% plotCorrResponseAccReactionTime: Calculate and plot the correlation between accuracy and reaction time
%
% Inputs:
%   - accuracyResults_updt: Table containing accuracy data
%   - T_patn_MEGfMRI: Table containing reaction time data

% Extract relevant data
animalAcc = accuracyResults_updt.animalAcc;
falsefontAcc = accuracyResults_updt.symbolAcc;
animalRT = T_patn_MEGfMRI.animalRT;
falsefontRT = T_patn_MEGfMRI.symbolRT;

% Calculate correlation for animal condition
[correlationAnimal, pValueAnimal] = corr(animalAcc, animalRT, 'Rows', 'complete');

% Calculate correlation for falsefont condition
[correlationFalsefont, pValueFalsefont] = corr(falsefontAcc, falsefontRT, 'Rows', 'complete');

% Display results
fprintf('Animal Corr = %.3f, p = %.3f\n', correlationAnimal, pValueAnimal);
fprintf('Falsefont Corr = %.3f, p = %.3f\n', correlationFalsefont, pValueFalsefont);

% Plot
figure('Color','w');
set(gcf, 'Position', [1000, 400, 450, 450]);

% ============== Animal ==============
subplot(2,1,1);
hold on;
scatter_size = 40;

scatter(animalAcc, animalRT, scatter_size, [0.2, 0.6, 0.8], 'filled');
scatter(animalAcc(discordSubs), animalRT(discordSubs), ...
        scatter_size, ...
        'r', ...
        'o', ...
        'LineWidth',1.5, ...
        'MarkerFaceColor','none');  % or 'filled' if you prefer

% Fit a line
pAnimal = polyfit(animalAcc, animalRT, 1);
xFitA = linspace(min(animalAcc), max(animalAcc), 100);
yFitA = polyval(pAnimal, xFitA);
plot(xFitA, yFitA, 'Color', [.6 .6 .6], 'LineWidth', 2);

title('Animal: Acc vs Reaction Time');
xlabel('Animal Accuracy');
ylabel('Reaction Time (s)');
grid on;
box off;
set(gca, 'color', 'none');
hold off;

% ============== Falsefont ==============
subplot(2,1,2);
hold on;
hScatter1 = scatter(falsefontAcc, falsefontRT, scatter_size, [0.2, 0.6, 0.8], 'filled', 'DisplayName','All subjects');
hScatter2 = scatter(falsefontAcc(discordSubs), falsefontRT(discordSubs), ...
        scatter_size, ...
        'r', ...
        'o', ...
        'LineWidth',1.5, ...
        'MarkerFaceColor','none',...
        'DisplayName','Discordant');  % or 'filled' if you prefer

legend([hScatter1, hScatter2], 'Location','best'); 

title('Falsefont: Acc vs Reaction Time');
xlabel('Falsefont Accuracy');
ylabel('Reaction Time (s)');
grid on;
box off;
set(gca, 'color', 'none');
hold off;

% figure
% % Fit a line
% pFF = polyfit(falsefontAcc, falsefontRT, 1);
% xFitFF = linspace(min(falsefontAcc), max(falsefontAcc), 100);
% yFitFF = polyval(pFF, xFitFF);
% plot(xFitFF, yFitFF, 'Color', [.6 .6 .6], 'LineWidth', 2);

end

