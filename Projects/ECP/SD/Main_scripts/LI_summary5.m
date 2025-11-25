% % -----------------------------------------------------------
% % Concordance bar charts per ROI:
% % Methods × (A-vs-B/A-vs-S, Fixed/Peak)
% % -----------------------------------------------------------
% 
% close all
% 
% methods   = {'Magnitude','Counting','Bootstrap'};
% 
% barLabels = { 'SD-vs-Base | Fixed interval' ...
%               'SD-vs-Base | Patient-Specific'  ...
%               'SD-vs-SYM | Fixed interval' ...
%               'SD-vs-SYM | Patient-Specific'  };
% 
% % Colour-blind-safe palette (4 columns)
% cmap = [ 0.96 0.60 0.15 ;   % orange
%          0.85 0.40 0.18 ;   % dark orange
%          0.40 0.76 0.65 ;   % teal
%          0.23 0.50 0.70 ];  % blue-teal
% 
% % --------- Data from the ROI table (rows = methods, cols = 4 bar types) -----
% % Angular
% M_ang = [ ...
%     69.0   66.2   81.9   80.5 ;   % Magnitude
%     69.0   64.79  81.9   83.3 ;   % Vertex counting
%     67.6   63.3   83.3   81.9 ];  % Bootstrap
% 
% % Frontal
% M_fro = [ ...
%     71.8   66.2   80.5   86.1 ;
%     71.8   66.2   83.3   80.5 ;
%     71.8   64.7   84.7   83.3 ];
% 
% % Temporal
% M_tmp = [ ...
%     63.3   66.2   84.7   84.7 ;
%     63.3   69.0   84.7   81.9 ;
%     63.3   63.3   83.3   86.1 ];
% 
% % Lateral
% M_lat = [ ...
%     71.8   73.2   88.8   88.8 ;
%     71.8   71.8   86.1   87.5 ;
%     67.6   70.4   87.5   90.2 ];
% 
% M_all   = {M_ang, M_fro, M_tmp, M_lat};
% roiNames = {'Angular','Frontal','Temporal','Lateral'};
% 
% figure('Units','pixels','Position',[100 100 400 1000])
% 
% for r = 1:4
%     subplot(4,1, r)
%     M = M_all{r};
% 
%     b = bar(M,'grouped'); hold on
%     for k = 1:numel(b)
%         b(k).FaceColor = cmap(k,:);
%         b(k).EdgeColor = 'none';
%     end
% 
%     ylim([60 100]); yticks(60:10:100); grid on
%     ax = gca; 
%     ax.GridLineStyle=':'; ax.YGrid='on';
%     ax.XTick = 1:numel(methods); 
%     ax.XTickLabel = methods;
%     ax.FontSize = 9; 
%     ax.Box = 'off';
%     set(gca,'color','none');
% 
%     title(roiNames{r},'FontSize',11,'FontWeight','normal')
% 
%     if r == 4
%         ylabel('LI concordance (%)','FontSize',10)
%     end
% 
%     % Value labels inside bars
% %     for i = 1:size(M,1)
% %         for j = 1:size(M,2)
% %             x = b(j).XData(i) + b(j).XOffset;
% %             y = M(i,j);
% %             text(x, y/2, sprintf('%.1f',y), ...
% %                  'HorizontalAlignment','center', ...
% %                  'VerticalAlignment','middle', ...
% %                  'Color','k','FontSize',7,'FontWeight','bold')
% %         end
% %     end
% 
%     % Put the legend only on the first subplot
%     if r == 4
%         lg = legend(barLabels,'Location','northoutside', ...
%                     'Orientation','horizontal','FontSize',8,'Box','off');
%     end
% end

data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim/';
cd(data_save_dir)

close all

M_all   = {M_ang, M_fro, M_tmp, M_lat};
roiNames = {'Angular','Frontal','Temporal','Lateral'};

% ===== Main ROI figure =====
hFig = figure('Units','pixels','Position',[100 100 400 1000]);

for r = 1:4
    subplot(4,1, r)
    M = M_all{r};

    b = bar(M,'grouped'); hold on
    for k = 1:numel(b)
        b(k).FaceColor = cmap(k,:);
        b(k).EdgeColor = 'none';
    end

    % Save bar handles from the FIRST subplot for legend later
    if r == 1
        bLegend = b;   % store handles
    end

    ylim([60 100]); yticks(60:10:100); grid on
    ax = gca; 
    ax.GridLineStyle=':'; ax.YGrid='on';
    ax.XTick = 1:numel(methods); 
    ax.XTickLabel = methods;
    ax.FontSize = 9; 
    ax.Box = 'off';
    set(gca,'color','none');

    title(roiNames{r},'FontSize',11,'FontWeight','normal')

    if r == 4
        ylabel('LI concordance (%)','FontSize',10)
    end
end

% ===== Separate legend figure (2x2 layout) =====
hLeg = figure('Units','pixels','Position',[100 100 500 100]);
axL = axes('Position',[0 0 1 1],'Visible','off');  % dummy axes

lg = legend(axL, bLegend, barLabels, ...
    'Orientation','horizontal', ...
    'NumColumns', 2, ...
    'Box','off', ...
    'FontSize',9);
axis off

% ===== Save both figures as SVG =====
main_svg_name = fullfile(data_save_dir, 'Fig_LI_concordance_ROIs.svg');
leg_svg_name  = fullfile(data_save_dir, 'Fig_LI_concordance_legend.svg');

print(hFig, main_svg_name, '-dsvg', '-painters');
print(hLeg, leg_svg_name,  '-dsvg', '-painters');

% Alternatively (older MATLAB):
% saveas(hFig, main_svg_name, 'svg');
% saveas(hLeg, leg_svg_name,  'svg');


%%
% Assume M_ang, M_fro, M_tmp, M_lat, methods already exist

% Stack into 3 x 4 x 4 (methods x barTypes x ROI)
M_all = cat(3, M_ang, M_fro, M_tmp, M_lat);

% -------- Overall stability (all ROIs + all 4 bar types) ----------
M_flat = reshape(M_all, 3, []);   % 3 x 16

mean_all = mean(M_flat, 2);
std_all  = std(M_flat, 0, 2);     % SD as stability measure

fprintf('\nOverall (all ROIs & contrasts):\n');
for m = 1:numel(methods)
    fprintf('%s: mean = %.2f, SD = %.2f\n', ...
        methods{m}, mean_all(m), std_all(m));
end

% -------- Stability for SD-vs-Baseline (cols 12) ----------
M_base = M_all(:, 1:2, :);             % 3 x 2 x 4
M_base_flat = reshape(M_base, 3, []);  % 3 x 8
mean_base = mean(M_base_flat, 2);
std_base  = std(M_base_flat, 0, 2);

fprintf('\nSD vs. Baseline only:\n');
for m = 1:numel(methods)
    fprintf('%s: mean = %.2f, SD = %.2f\n', ...
        methods{m}, mean_base(m), std_base(m));
end

% -------- Stability for SD-vs-SYM (cols 34) ----------
M_sym = M_all(:, 3:4, :);              % 3 x 2 x 4
M_sym_flat = reshape(M_sym, 3, []);    % 3 x 8
mean_sym = mean(M_sym_flat, 2);
std_sym  = std(M_sym_flat, 0, 2);

fprintf('\nSD vs. SYM only:\n');
for m = 1:numel(methods)
    fprintf('%s: mean = %.2f, SD = %.2f\n', ...
        methods{m}, mean_sym(m), std_sym(m));
end
