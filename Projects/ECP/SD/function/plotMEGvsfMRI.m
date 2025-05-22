function demingRes = plotMEGvsfMRI(optMEG_LI, fMRI_LI, discordSubs, subIDs, lambda)
% PLOTMEGVSFMRI  Scatter of MEG vs fMRI laterality indices with Deming fit.
%
% demingRes = plotMEGvsfMRI(optMEG_LI, fMRI_LI, discordSubs, subIDs, lambda)
%
% INPUTS
%   optMEG_LI   [N×1] MEG LIs
%   fMRI_LI     [N×1] fMRI LIs
%   discordSubs vector of indices for discordant subjects (can be empty [])
%   subIDs      [N×1] numeric OR cell array of char labels
%   lambda      scalar s²_y / s²_x (default = 1 ? equal error variance)
%
% OUTPUT
%   demingRes   struct with Deming fit results:
%                 .b0, .b1, .b0_CI, .b1_CI, .res_SE, .sigma2x, .lambda
%
% The plot shows:
%    All subjects (blue), discordant (red)
%    y = x reference line
%    Deming line ±95 % CI (magenta solid / dashed)
%    Pearson r

% ----------------------- sanity checks ----------------------------------
narginchk(4,5)
if numel(optMEG_LI) ~= numel(fMRI_LI) || numel(optMEG_LI) ~= numel(subIDs)
    error('optMEG_LI, fMRI_LI, and subIDs must have the same length.')
end
if nargin < 5 || isempty(lambda);  lambda = 1;  end
if isempty(discordSubs);           discordSubs = []; end
% ------------------------------------------------------------------------

%% ---- figure & scatter -------------------------------------------------
figure('Color','w','Name','MEG vs fMRI LI','Position',[220 180 720 520]);
scatter(optMEG_LI, fMRI_LI, 50, 'o', ...
        'MarkerFaceColor',[0.2 0.6 0.8], 'MarkerEdgeColor','none', ...
        'DisplayName','All subjects');  hold on; grid on

if ~isempty(discordSubs)
    scatter(optMEG_LI(discordSubs), fMRI_LI(discordSubs), 50, ...
            'MarkerFaceColor','r','MarkerEdgeColor','none', ...
            'DisplayName','Gross discordant');
end

xlabel('Right <                MEG LI                > Left')
ylabel('Right <                fMRI LI               > Left')
xlim([-100 100]); ylim([-100 100])
set(gca,'XTick',-100:50:100,'YTick',-100:50:100,'FontSize',10)
axis square

% subject labels ---------------------------------------------------------
for k = 1:numel(subIDs)
    text(optMEG_LI(k), fMRI_LI(k), num2str(subIDs(k)), ...
        'VerticalAlign','bottom','HorizontalAlign','left','FontSize',8);
end

% y = x reference line ---------------------------------------------------
plot(xlim, ylim, 'k--', 'LineWidth',1, 'DisplayName','y = x')

% Pearson correlation ----------------------------------------------------
[r,p]   = corr(optMEG_LI(:), fMRI_LI(:), 'Rows','complete');
text(-95, 90, sprintf('r = %.2f (p = %.3g)', r, p), 'Color','b', 'FontSize',10)

%% ---------- Deming regression line ------------------------------------
[b, sigma2_x, ~, ~, stats] = deming(optMEG_LI, fMRI_LI, lambda);

xFit  = linspace(min(optMEG_LI), max(optMEG_LI), 100);
yFit  = b(1) + b(2)*xFit;
% plot(xFit, yFit, 'm-', 'LineWidth',1.8, 'DisplayName','Deming fit')

% if isfield(stats,'b_ci') && ~isempty(stats.b_ci)  % CI bands if available
%     yLo = stats.b_ci(1,1) + stats.b_ci(2,1)*xFit;
%     yHi = stats.b_ci(1,2) + stats.b_ci(2,2)*xFit;
%     plot(xFit, yLo, 'm--', xFit, yHi, 'm--', 'HandleVisibility','off')
% end

%% --------- package results for manuscript ------------------------------
demingRes = struct( ...
    'b0',       b(1), ...
    'b1',       b(2), ...
    'b0_CI',    stats.b_ci(1,:), ...
    'b1_CI',    stats.b_ci(2,:), ...
    'res_SE',   stats.s_e, ...
    'sigma2x',  sigma2_x, ...
    'lambda',   lambda);

legend('Location','bestoutside')


%% ---------- compact Deming summary on the plot -------------------------
% Format numbers
% fmtNum   = @(x) sprintf('%.2f', x);          % two decimals
% fmtIntvl = @(ci) sprintf('[%s , %s]', fmtNum(ci(1)), fmtNum(ci(2)));
% 
% sumLines = { ...
%     sprintf('\\bfDeming fit (\\lambda = %g)', demingRes.lambda)   ...
%     , sprintf('Intercept  b_0 = %s  %s', fmtNum(demingRes.b0), ...
%               fmtIntvl(demingRes.b0_CI))                         ...
%     , sprintf('Slope      b_1 = %s  %s', fmtNum(demingRes.b1), ...
%               fmtIntvl(demingRes.b1_CI))                         ...
%     , sprintf('Residual  s_e = %.1f LI units', demingRes.res_SE) ...
%     };
% 
% % Place as a textbox in normalised axes units
% annotation('textbox', [0.62 0.02 0.35 0.18], ...  % [x y w h] in figure
%            'String',  sumLines, ...
%            'Interpreter','tex', ...
%            'FontSize', 9, ...
%            'EdgeColor', 'none', ...
%            'BackgroundColor', [1 1 1 0.75]);      % white, slightly transparent



hold off



end


% function plotMEGvsfMRI(optMEG_LI, fMRI_LI, discordSubs, subIDs)
% % PLOTMEGVSFMRI  Plots a scatter of MEG LI vs. fMRI LI, highlights discordant subjects in red,
% %                and includes subject labels plus correlation/regression lines.
% %
% %   plotMEGvsfMRI(optMEG_LI, fMRI_LI, discordSubs, subIDs)
% %
% %   INPUTS:
% %       optMEG_LI     - Nx1 numeric array of MEG LI values
% %       fMRI_LI       - Nx1 numeric array of fMRI LI values (same length as optMEG_LI)
% %       discordSubs   - Index vector of "discordant" subjects (e.g., [4,5,10,...])
% %       subIDs        - Nx1 numeric or cell array for subject labels (same length as optMEG_LI).
% %                       If you don't have custom labels, pass 1:N or something similar.
% %
% %   The function:
% %       1) Creates a figure, plots all points in a custom color.
% %       2) Overlays discordant samples in red.
% %       3) Labels each point with subIDs.
% %       4) Sets axis labels to indicate right (<0) vs. left (>0).
% %       5) Sets axis limits and ticks to [-100,100] by 50 steps (customizable).
% %       6) Plots y=x line, computes correlation, and plots a simple linear fit.
% %       7) Adds a legend and styling adjustments.
% 
% if nargin < 4
%     error('All 4 inputs are required: optMEG_LI, fMRI_LI, discordSubs, subIDs.');
% end
% if length(optMEG_LI) ~= length(fMRI_LI) || length(optMEG_LI) ~= length(subIDs)
%     error('Input arrays must have the same length: optMEG_LI, fMRI_LI, subIDs.');
% end
% 
% % Create figure
% figure('Color','w','Name','MEG vs fMRI LI','Position',[200,200,700,500]);
% 
% % 1) Plot all subjects in one color
% scatter(optMEG_LI, fMRI_LI, 50, 'o',...
%     'MarkerFaceColor',[0.2, 0.6, 0.8],...
%     'MarkerEdgeColor','none',...
%     'DisplayName','All Subjects');
% hold on; grid on;
% 
% % 2) Overlay discordant subjects in red
% scatter(optMEG_LI(discordSubs), fMRI_LI(discordSubs), 50,...
%     'MarkerFaceColor','r',...
%     'MarkerEdgeColor','none',...
%     'DisplayName','Discordant');
% 
% % 3) Axis labels indicating negative=Right, positive=Left
% xlabel('Right <                MEG LI                > Left');
% ylabel('Right <                fMRI LI               > Left');
% 
% % 4) Set axis limits to [-100, 100] and matching ticks
% xlim([-100, 100]); ylim([-100, 100]);
% set(gca, 'XTick', -100:50:100, 'YTick', -100:50:100);
% 
% % 5) Label each subject (optional)
% for s = 1:length(optMEG_LI)
%     text(optMEG_LI(s), fMRI_LI(s), num2str(subIDs(s)), ...
%         'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
%         'FontSize',8, 'Color','k');
% end
% 
% % 6) Plot y=x reference line
% xLimits = xlim;
% plot(xLimits, xLimits, 'k--', 'LineWidth',1, 'DisplayName','y = x');
% 
% % 7) Compute correlation & show in plot
% [r, p] = corr(optMEG_LI, fMRI_LI, 'Rows','complete');
% corrStr = sprintf('r = %.2f (p=%.3g)', r, p);
% text(xLimits(1)+5, xLimits(2)-5, corrStr, 'FontSize',10, 'Color','b');
% 
% % 8) Simple linear fit
% % pCoeffs = polyfit(optMEG_LI, fMRI_LI, 1);  % slope/intercept
% % xVals   = linspace(min(optMEG_LI), max(optMEG_LI), 100);
% % yFit    = polyval(pCoeffs, xVals);
% % plot(xVals, yFit, 'b-', 'LineWidth',1.5, 'DisplayName','Corr line');
% 
% %% ----- Replace the LS line with a Deming fit ---------------------------
% % slope/intercept (b = [b0; b1])
% lambda = 1;
% [b,sigma2_x,~,~,stats] = deming(optMEG_LI, fMRI_LI, lambda);   % <-- your deming.m
% % Draw Deming line
% xVals = linspace(min(optMEG_LI), max(optMEG_LI), 100);
% yDem  = b(1) + b(2)*xVals;
% plot(xVals, yDem, 'm-', 'LineWidth',1.8, 'DisplayName','Deming fit');
% 
% % Optional 95 % CI (needs Statistics Toolbox for tinv)
% if isfield(stats,'b_ci') && ~isempty(stats.b_ci)
%     b_lo = stats.b_ci(:,1);   b_hi = stats.b_ci(:,2);
%     y_lo = b_lo(1) + b_lo(2)*xVals;
%     y_hi = b_hi(1) + b_hi(2)*xVals;
%     plot(xVals, y_lo,'m--', xVals, y_hi,'m--','HandleVisibility','off');
% end
% % -----------------------------------------------------------------------
% 
% % Convenience struct for table or save
% demingRes.b0      = b(1);
% demingRes.b1      = b(2);
% demingRes.b0_CI   = stats.b_ci(1,:);
% demingRes.b1_CI   = stats.b_ci(2,:);
% demingRes.res_SE  = stats.s_e;
% demingRes.sigma2x = sigma2_x;
% demingRes.lambda  = 1;
% 
% 
% % 9) Legend & styling
% legend('Location','bestoutside');
% axis tight;
% set(gca,'FontName','Helvetica','FontSize',10);
% axis square;
% 
% hold off;
% end
