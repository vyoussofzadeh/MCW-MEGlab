
%% Resources
% 1) https://brendanhasz.github.io/2019/07/03/matlab-uncertainty-viz.html

%% Matlab figure settings
figure,
bar((DL_runs.math + DL_runs.str)/2), L = length(DL_runs.math);
set(gca,'Xtick', 1:L,'XtickLabel',sub_runs_name);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   100   1500   300]);
set(gca,'color','none');
title([tag, ', DL, mean (Str, Math) '])
axis square
set(lgnd,'color','none');

%% Overlay text on bargraph
text(1:length(Y),Y,num2str(Y'),'vert','bottom','horiz','center');
text(1:length(Y),Y,num2str(Y','%4.2f'),'vert','bottom','horiz','center');

%% Best fit to a figure,
axis tight

%% Legend
legend('hide')

%% Axis labelling
xlabel(rois(idx(1)))
ylabel([tag,'-',ttag])

%% Corrolation
[R,P] = corrcoef(X(:,idx(1)),y1');

%% Regression
tbl = table(X(:,idx(1)),y1');
mdl = fitlm(tbl,'linear');
plot(mdl)

%% some helpful links
% https://anneurai.net/2016/06/13/prettier-plots-in-matlab/

%% Distinguishable_colors
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External');
colr = distinguishable_colors(nScouts);

%% Shaded area error bar plot
% link: 'https://brendanhasz.github.io/2019/07/03/matlab-uncertainty-viz.html'

N = 60;
X = linspace(-2, 6, N)';
Y1 = exp(-X)-exp(-2*X);
Y1(X<0) = 0;
Y1 = Y1 + 0.01*randn(N, 1);
E1 = 0.02+0.1*rand(N, 1);
Y2 = exp(-X+1)-exp(-2*(X-1));
Y2(Y2<0) = 0;
Y2 = Y2 + 0.01*randn(N, 1);
E2 = 0.03+0.05*rand(N, 1);

% Plot shaded, semitransparent error bounds
blue = [0.35 0.7 0.9]; orange = [0.9,0.6,0];

figure
fill([X; flipud(X)], [Y1+E1; flipud(Y1-E1)], blue, ...
    'EdgeColor', 'none', 'facealpha', 0.3)
hold on
plot(X, Y1, 'Color', blue, 'LineWidth', 2)
fill([X; flipud(X)], [Y2+E2; flipud(Y2-E2)], orange, ...
    'EdgeColor', 'none', 'facealpha', 0.3)
plot(X, Y2, 'Color', orange, 'LineWidth', 2)

% example 2
figure,
y = rand(1,10); % your mean vector;
x = 1:numel(y);
std_dev = 1;
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'g');
hold on;
plot(x, y, 'r', 'LineWidth', 2);

%% Legned stop
legend('AutoUpdate', 'off')

h2 = fill(x2, inBetween, colr(j,:), 'EdgeColor', 'none', 'facealpha', 0.1);
h2.Annotation.LegendInformation.IconDisplayStyle = 'off'; % make the legend for step plot off

%%