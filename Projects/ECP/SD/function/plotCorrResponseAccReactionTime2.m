function plotCorrResponseAccReactionTime2(acc, rt)
% plotCorrResponseAccReactionTime: Calculate and plot the correlation between accuracy and reaction time
%
% Inputs:
%   - accuracyResults_updt: Table containing accuracy data
%   - T_patn_MEGfMRI: Table containing reaction time data

% Convert numerical subject IDs in accuracyResults_updt to string format with 'EC' prefix
% accuracyResults_updt.SubjectStr = strcat('EC', arrayfun(@num2str, accuracyResults_updt.Subject, 'UniformOutput', false));

animalAcc = acc.anim;
falsefontAcc = acc.falsefont;

animalRT = rt.anim;
falsefontRT = rt.falsefont;

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
scatter(falsefontAcc, falsefontRT, scatter_size, [0.2, 0.6, 0.8], 'filled');

% Fit a line
pFF = polyfit(falsefontAcc, falsefontRT, 1);
xFitFF = linspace(min(falsefontAcc), max(falsefontAcc), 100);
yFitFF = polyval(pFF, xFitFF);
plot(xFitFF, yFitFF, 'Color', [.6 .6 .6], 'LineWidth', 2);

title('Falsefont: Acc vs Reaction Time');
xlabel('Falsefont Accuracy');
ylabel('Reaction Time (s)');
grid on;
box off;
set(gca, 'color', 'none');
hold off;

end

