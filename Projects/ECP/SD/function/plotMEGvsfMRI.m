function plotMEGvsfMRI(optMEG_LI, fMRI_LI, discordSubs, subIDs)
% PLOTMEGVSFMRI  Plots a scatter of MEG LI vs. fMRI LI, highlights discordant subjects in red,
%                and includes subject labels plus correlation/regression lines.
%
%   plotMEGvsfMRI(optMEG_LI, fMRI_LI, discordSubs, subIDs)
%
%   INPUTS:
%       optMEG_LI     - Nx1 numeric array of MEG LI values
%       fMRI_LI       - Nx1 numeric array of fMRI LI values (same length as optMEG_LI)
%       discordSubs   - Index vector of "discordant" subjects (e.g., [4,5,10,...])
%       subIDs        - Nx1 numeric or cell array for subject labels (same length as optMEG_LI).
%                       If you don't have custom labels, pass 1:N or something similar.
%
%   The function:
%       1) Creates a figure, plots all points in a custom color.
%       2) Overlays discordant samples in red.
%       3) Labels each point with subIDs.
%       4) Sets axis labels to indicate right (<0) vs. left (>0).
%       5) Sets axis limits and ticks to [-100,100] by 50 steps (customizable).
%       6) Plots y=x line, computes correlation, and plots a simple linear fit.
%       7) Adds a legend and styling adjustments.

if nargin < 4
    error('All 4 inputs are required: optMEG_LI, fMRI_LI, discordSubs, subIDs.');
end
if length(optMEG_LI) ~= length(fMRI_LI) || length(optMEG_LI) ~= length(subIDs)
    error('Input arrays must have the same length: optMEG_LI, fMRI_LI, subIDs.');
end

% Create figure
figure('Color','w','Name','MEG vs fMRI LI','Position',[200,200,700,500]);

% 1) Plot all subjects in one color
scatter(optMEG_LI, fMRI_LI, 50, 'o',...
    'MarkerFaceColor',[0.2, 0.6, 0.8],...
    'MarkerEdgeColor','none',...
    'DisplayName','All Subjects');
hold on; grid on;

% 2) Overlay discordant subjects in red
scatter(optMEG_LI(discordSubs), fMRI_LI(discordSubs), 50,...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','none',...
    'DisplayName','Discordant');

% 3) Axis labels indicating negative=Right, positive=Left
xlabel('Right <                MEG LI                > Left');
ylabel('Right <                fMRI LI               > Left');

% 4) Set axis limits to [-100, 100] and matching ticks
xlim([-100, 100]); ylim([-100, 100]);
set(gca, 'XTick', -100:50:100, 'YTick', -100:50:100);

% 5) Label each subject (optional)
for s = 1:length(optMEG_LI)
    text(optMEG_LI(s), fMRI_LI(s), num2str(subIDs(s)), ...
        'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
        'FontSize',8, 'Color','k');
end

% 6) Plot y=x reference line
xLimits = xlim;
plot(xLimits, xLimits, 'k--', 'LineWidth',1, 'DisplayName','y = x');

% 7) Compute correlation & show in plot
[r, p] = corr(optMEG_LI, fMRI_LI, 'Rows','complete');
corrStr = sprintf('r = %.2f (p=%.3g)', r, p);
text(xLimits(1)+5, xLimits(2)-5, corrStr, 'FontSize',10, 'Color','b');

% 8) Simple linear fit
pCoeffs = polyfit(optMEG_LI, fMRI_LI, 1);  % slope/intercept
xVals   = linspace(min(optMEG_LI), max(optMEG_LI), 100);
yFit    = polyval(pCoeffs, xVals);
plot(xVals, yFit, 'b-', 'LineWidth',1.5, 'DisplayName','Corr line');

% 9) Legend & styling
legend('Location','bestoutside');
axis tight;
set(gca,'FontName','Helvetica','FontSize',10);
axis square;

hold off;
end
