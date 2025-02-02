function ecpfunc_blandAltmanPlot(megLI, fmriLI, subjectLabels)
% ECPFUNC_BLANDALTMANPLOT_OUTLIERSONLY
%   Creates a BlandAltman Plot for MEG LI vs. fMRI LI, labeling only
%   those data points (subjects) whose difference is outside ±1.96 SD,
%   and showing the mean diff / limits in a separate textbox instead of on
%   the lines themselves (to avoid obscuring points).
%
% INPUTS:
%   megLI         : [N x 1] or [1 x N], MEG lateralization indices
%   fmriLI        : [N x 1] or [1 x N], fMRI lateralization indices
%   subjectLabels : cell array of strings, length N (e.g. {'S1','S2','S3',...})

    % Ensure column vectors
    megLI  = megLI(:);
    fmriLI = fmriLI(:);
    nSubj  = length(megLI);

    if nargin < 3 || isempty(subjectLabels)
        % If no labels provided, make dummy labels or none
        subjectLabels = arrayfun(@(x) sprintf('S%d', x), 1:nSubj, 'UniformOutput', false);
    elseif length(subjectLabels) ~= nSubj
        error('subjectLabels must match the length of megLI/fmriLI.');
    end

    % 1) Compute BlandAltman metrics
    meanVals = (megLI + fmriLI) / 2;
    diffVals = megLI - fmriLI;
    meanDiff = mean(diffVals);
    sdDiff   = std(diffVals);
    loaUpper = meanDiff + 1.96 * sdDiff;
    loaLower = meanDiff - 1.96 * sdDiff;

    % 2) Create figure & scatter
    figure('Color','w'); 
    scatter(meanVals, diffVals, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on; grid on; box on;

    % 3) Horizontal lines (no text labels on lines, to avoid obscuring data)
    yline(meanDiff, 'r--', 'LineWidth', 2);
    yline(loaUpper, 'k--', 'LineWidth', 1.5);
    yline(loaLower, 'k--', 'LineWidth', 1.5);

    xlabel('Mean of MEG LI and fMRI LI');
    ylabel('Difference (MEG LI - fMRI LI)');
    title('Bland-Altman Plot');

    % 4) Show stats in a textbox (outside the plot region, e.g., upper left corner)
    statsString = { ...
       sprintf('Mean Diff = %.3f', meanDiff), ...
       sprintf('+1.96 SD = %.3f', loaUpper), ...
       sprintf('-1.96 SD = %.3f', loaLower) ...
    };
    annotation('textbox', [0.15 0.82 0.1 0.1], ...  % [x y width height] in normalized coords
        'String', statsString, ...
        'FitBoxToText', 'on', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'white', ...
        'FontSize', 9);

    % 5) Label ONLY outliers (points outside ±1.96 SD)
    for i = 1:nSubj
        if diffVals(i) > loaUpper || diffVals(i) < loaLower
            text(meanVals(i), diffVals(i), subjectLabels{i}, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment',   'bottom', ...
                'FontSize', 8, ...
                'Color', 'blue');
        end
    end

end
