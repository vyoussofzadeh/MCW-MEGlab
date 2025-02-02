function plotConcMaxByROI(summaryTable, baseColors)
% PLOTCONCMAXBYROI Plots "Concordance" maximum (Max_Value) by ROI using a
%                  color-darkening scheme for higher values.
%
%   plotConcMaxByROI(summaryTable) plots bar charts for each ROI in
%   summaryTable, using a default set of base colors.
%
%   plotConcMaxByROI(summaryTable, baseColors) allows you to specify
%   a custom Nx3 matrix of RGB base colors. If there are more ROIs
%   than rows in baseColors, the colors are reused in a round-robin fashion.
%
%   REQUIRED COLUMNS in summaryTable (type: table):
%       - ROI (categorical or string)
%       - Metric_Type (categorical or string) with 'Concordance' entries
%       - LI_Method (categorical or string)
%       - Max_Value (numeric)
%
%   The function automatically:
%       1. Segregates ROI data and filters rows for Metric_Type == 'Concordance'
%       2. Normalizes Max_Value in [0,1] to compute an inverted darkening factor
%       3. Produces one subplot per ROI, with darker bars for higher values.
%   Bars are narrower, and the figure is sized to be smaller by default.
%
%   Example:
%       % Basic Usage
%       plotConcMaxByROI(mySummaryTable);
%
%       % Custom Base Colors (4 example colors)
%       customColors = [
%           0.96, 0.49, 0.00;  % Orange
%           0.22, 0.56, 0.24;  % Green
%           0.69, 0.71, 0.17;  % Yellow
%           0.48, 0.12, 0.66   % Purple
%       ];
%       plotConcMaxByROI(mySummaryTable, customColors);

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
        error('The summaryTable is missing required variables: %s', strjoin(missingVars, ', '));
    end

    %% 2. Prepare Data & Figure
%     uniqueROIs = unique(summaryTable.ROI);
    uniqueROIs = {'Ang' 'Front' 'Temp' 'Lat'};
    nROIs      = numel(uniqueROIs);

    % Create the main figure (white background), smaller size
    figure('Color','w','Units','normalized','Position',[0.3 0.2 0.10 0.5]);
    
    % If sgtitle is available (R2018b+), use it:
    if exist('sgtitle','file') || ~isempty(which('sgtitle'))
        sgtitle('Concordance Max by ROI','FontSize',12,'FontWeight','bold');
    else
        % For older MATLAB versions, fallback to suptitle if available
        if exist('suptitle','file')
            suptitle('Concordance Max by ROI');
        else
            set(gcf,'Name','Concordance Max by ROI');
        end
    end

    %% 3. Loop Over ROIs & Create One Subplot per ROI
    for i = 1:nROIs
        % Create a subplot for this ROI
        subplot(nROIs, 1, i);
        
        roi = uniqueROIs{i};
        
        % Pick (or cycle) the base color for this ROI
        colorIndex = mod(i-1, size(baseColors, 1)) + 1;
        baseColor  = baseColors(colorIndex, :);

        % Extract data for the current ROI, filter for 'Concordance'
        roiData  = summaryTable(strcmp(summaryTable.ROI, roi), :);
        concData = roiData(strcmp(roiData.Metric_Type, 'Concordance'), :);

        % Handle empty data
        if isempty(concData)
            title(['ROI: ', char(roi), ' (No Concordance Data)'], 'FontWeight','bold');
            ylim([0,1]);
            continue;
        end

        %% 4. Normalize Data [0,1] & Invert Scale for Darker = Higher
        maxVal = max(concData.Max_Value);
        minVal = min(concData.Max_Value);
        if maxVal == minVal
            colorScaling = zeros(size(concData.Max_Value));
        else
            colorScaling = (concData.Max_Value - minVal) ./ (maxVal - minVal);
        end

        % For higher values => smaller multiplier => darker color
        darkeningFactor = 1.0 - 0.5 * colorScaling;  % [1..0.5]

        %% 5. Plot Bars for each Method
        methods = unique(concData.LI_Method, 'stable');
        hold on;
        for j = 1:length(methods)
            method = methods{j};
            value  = concData.Max_Value(strcmp(concData.LI_Method, method));

            % Determine shading based on normalized value
            thisDarkening = darkeningFactor(strcmp(concData.LI_Method, method));
            scaledColor   = baseColor .* thisDarkening;

            % Plot bar with a narrower width (e.g., 0.2)
            bar(categorical({method}), value, ...
                'BarWidth', 0.2, ...
                'FaceColor', scaledColor, ...
                'EdgeColor', 'none');

            % Add numeric label above each bar
            text(j, value+2, sprintf('%.2f', value), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',9,'FontWeight','bold');
        end
        hold off;

        %% 6. Style Each Subplot
        title(['ROI: ', char(roi)], 'FontWeight','bold','FontSize',10);
        ylabel('Concordance (%)','FontWeight','bold','FontSize',9);

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
        set(gca,'YLim',[0, max(yLimits(2), 105)]);
    end

    %% 7. (Optional) Adjust Spacing
    % If you want to reduce or increase spacing between subplots, you can
    % manipulate 'Position' of each axis or 'PaperPositionMode'. For example:
%     for sp = 1:nROIs
%         ax = subplot(nROIs,1,sp);
%         pos = get(ax,'Position');
%         pos(4) = pos(4)*1.05;  % slightly taller
%         set(ax,'Position',pos);
%     end
end
