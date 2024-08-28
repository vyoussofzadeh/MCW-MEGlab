function plotCorrResponseAccReactionTime(accuracyResults_updt, T_patn_MEGfMRI)
    % plotCorrResponseAccReactionTime: Calculate and plot the correlation between accuracy and reaction time
    %
    % Inputs:
    %   - accuracyResults_updt: Table containing accuracy data
    %   - T_patn_MEGfMRI: Table containing reaction time data

    % Convert numerical subject IDs in accuracyResults_updt to string format with 'EC' prefix
    accuracyResults_updt.SubjectStr = strcat('EC', arrayfun(@num2str, accuracyResults_updt.Subject, 'UniformOutput', false));

    % Find the intersection between subject IDs
    [~, idxAccuracy, idxRT] = intersect(accuracyResults_updt.SubjectStr, T_patn_MEGfMRI.Sub_ID);

    % Extract relevant data
    animalAcc = accuracyResults_updt.Animal_ACC(idxAccuracy);
    falsefontAcc = accuracyResults_updt.Falsefont_ACC(idxAccuracy);
    animalRT = T_patn_MEGfMRI.Animal(idxRT);
    falsefontRT = T_patn_MEGfMRI.Symbol(idxRT);

    % Calculate correlation for animal condition
    [correlationAnimal, pValueAnimal] = corr(animalAcc, animalRT, 'Rows', 'complete');

    % Calculate correlation for falsefont condition
    [correlationFalsefont, pValueFalsefont] = corr(falsefontAcc, falsefontRT, 'Rows', 'complete');

    % Display results
    disp(['Correlation between Animal Accuracy and Reaction Time: ', num2str(correlationAnimal)]);
    disp(['P-value for Animal Accuracy and Reaction Time: ', num2str(pValueAnimal)]);

    disp(['Correlation between Falsefont Accuracy and Reaction Time: ', num2str(correlationFalsefont)]);
    disp(['P-value for Falsefont Accuracy and Reaction Time: ', num2str(pValueFalsefont)]);

    % Plot Correlations

    % Plot for Animal condition
    figure;
    subplot(2, 1, 1)
    scatter(animalAcc, animalRT, 'filled');
    title('Correlation (Animal Acc, Reaction Time)');
    xlabel('Animal Acc');
    ylabel('Reaction Time (s)');
    grid on;
    set(gca, 'color', 'none');

    % Plot for Falsefont condition
    subplot(2, 1, 2)
    scatter(falsefontAcc, falsefontRT, 'filled');
    title('Correlation (Falsefont Acc, Reaction Time)');
    xlabel('Falsefont Acc');
    ylabel('Reaction Time (s)');
    grid on;
    box off;
    set(gcf, 'Position', [1000, 100, 400, 400]);
    set(gca, 'color', 'none');
end
