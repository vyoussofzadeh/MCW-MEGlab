function plotCorrMaxByROI(summaryTable, baseColors)
% PLOTCORRMAXBYROI Plots correlation maximum (Max_Value) by ROI using a
%                  color-darkening scheme for higher values, with numeric
%                  x-axis and optionally adjustable bar spacing.
%
%   plotCorrMaxByROI(summaryTable) plots bar charts for each ROI in
%   summaryTable, using a default set of base colors.
%
%   plotCorrMaxByROI(summaryTable, baseColors) allows you to specify
%   a custom Nx3 matrix of RGB base colors. If there are more ROIs
%   than rows in baseColors, the colors are reused in a round-robin fashion.
%
%   REQUIRED COLUMNS in summaryTable (type: table):
%       - ROI (categorical or string)
%       - Metric_Type (categorical or string) with 'Correlation' entries
%       - LI_Method (categorical or string)
%       - Max_Value (numeric)
%
%   This function:
%       1. Segregates ROI data and filters rows for Metric_Type == 'Correlation'.
%       2. Normalizes Max_Value in [0,1] to compute an inverted darkening factor 
%          (higher => darker).
%       3. Produces one subplot per ROI, placing bars at integer x-values
%          (1, 2, 3, ...) for more precise control over spacing.
%       4. Darker bars indicate higher Max_Value within each ROI.
%
%   Example:
%       % Basic Usage
%       plotCorrMaxByROI(mySummaryTable);
%
%       % Custom Base Colors (4 example colors)
%       customColors = [
%           0.96, 0.49, 0.00;  % Orange
%           0.22, 0.56, 0.24;  % Green
%           0.69, 0.71, 0.17;  % Yellow
%           0.48, 0.12, 0.66   % Purple
%       ];
%       plotCorrMaxByROI(mySummaryTable, customColors);

    %% 1. Validate & Set Defaults
    if nargin < 2 || isempty(baseColors)
        % Default baseColors if none provided
        baseColors = [
            0.96, 0.49, 0.00;  % Orange
            0.22, 0.56, 0.24;  % Green
            0.69, 0.71, 0.17;  % Yellow
            0.48, 0.12, 0.66   % Purple
        ];
    end

    if ~istable(summaryTable)
        error('Input summaryTable must be a MATLAB table.');
    end

    requiredVars = {'ROI','Metric_Type','LI_Method','Max_Value'};
    missingVars  = setdiff(requiredVars, summaryTable.Properties.VariableNames);
    if ~isempty(missingVars)
        error('The summaryTable is missing required variables: %s', ...
            strjoin(missingVars, ', '));
    end

    %% 2. Prepare Data & Figure
    uniqueROIs = unique(summaryTable.ROI);
    nROIs      = numel(uniqueROIs);

    % Create the main figure
    figure('Color','w','Units','normalized','Position',[0.3 0.2 0.10 0.5]);

    % If sgtitle is available (R2018b+), use it:
    if exist('sgtitle','file') || ~isempty(which('sgtitle'))
        sgtitle('Corr Max by ROI','FontSize',12,'FontWeight','bold');
    else
        % Fallback for older MATLAB versions
        if exist('suptitle','file')
            suptitle('Corr Max by ROI');
        else
            set(gcf,'Name','Corr Max by ROI');
        end
    end

    %% 3. Loop Over ROIs & Create One Subplot per ROI
    for i = 1:nROIs
        subplot(nROIs, 1, i);
        roi = uniqueROIs{i};

        % Pick (or cycle) the base color for this ROI
        colorIndex = mod(i-1, size(baseColors, 1)) + 1;
        baseColor  = baseColors(colorIndex, :);

        % Extract data for the current ROI, filter for 'Correlation'
        roiData  = summaryTable(strcmp(summaryTable.ROI, roi), :);
        corrData = roiData(strcmp(roiData.Metric_Type, 'Correlation'), :);

        % Handle empty data
        if isempty(corrData)
            title(['ROI: ', char(roi), ' (No Correlation Data)'], 'FontWeight','bold');
            ylim([0, 1]);
            continue;
        end

        %% 4. Normalize Data [0,1] & Invert Scale for Darker = Higher
        maxVal = max(corrData.Max_Value);
        minVal = min(corrData.Max_Value);

        if maxVal == minVal
            % All values the same => skip dynamic scaling
            colorScaling = zeros(size(corrData.Max_Value));
        else
            colorScaling = (corrData.Max_Value - minVal) ./ (maxVal - minVal);
        end

        % For higher values => smaller multiplier => darker color
        darkeningFactor = 1.0 - 0.5 * colorScaling;  % range: [1..0.5]
        
        
        %% 5. Plot Bars for each Method using numeric x-axis
        methods = unique(corrData.LI_Method, 'stable');
        hold on;
        for j = 1:length(methods)
            method = methods{j};
            value  = corrData.Max_Value(strcmp(corrData.LI_Method, method));
            
            % Determine shading based on normalized value
            thisDarkening = darkeningFactor(strcmp(corrData.LI_Method, method));
            scaledColor   = baseColor .* thisDarkening;
            
            % Plot each bar at x = j, with adjustable BarWidth for spacing
            barWidth = 0.2;  % Adjust this between 0-1 to control horizontal spacing
            bar(categorical({method}), value, 'BarWidth', barWidth, 'FaceColor', scaledColor, 'EdgeColor','none');
            
            % Numeric label above each bar
            text(j, value + 0.02*maxVal, sprintf('%.2f', value), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',9,'FontWeight','bold');
        end
        hold off;

        %% 6. Style the Subplot
        title(['ROI: ', char(roi)], 'FontWeight','bold','FontSize',10);
        ylabel('Correlation','FontWeight','bold','FontSize',9);

        set(gca, ...
            'Box','off', ...
            'FontName','Helvetica', ...
            'FontSize',9, ...
            'LineWidth',1.2, ...
            'YGrid','on', ...
            'GridColor',[0.8, 0.8, 0.8], ...
            'GridAlpha',0.5, ...
            'Color','none');

%         axis tight;        
        yLimits = get(gca,'YLim');
        % Extend the top limit to avoid label overlap
        set(gca,'YLim',[0, max(yLimits(2), 1)]);
    end
end

