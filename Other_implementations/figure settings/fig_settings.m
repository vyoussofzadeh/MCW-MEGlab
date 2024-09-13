
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

%% Legend position
p1 = plot();
lgd = legend([p1, p2, p3], 'Location', 'best','Orientation', 'horizontal'); % Add a legend
lgdPos = lgd.Position; % Get current position
lgdPos(2) = lgdPos(2) - 0.09; % Move legend down
lgdPos(1) = lgdPos(1) - 0.05; % Move legend down
lgd.Position = lgdPos;

%% Change legend position
lgd = legend(resultsTable.Method, 'Location', 'southoutside', 'NumColumns', 2, 'Orientation', 'horizontal');
lgdPos = lgd.Position; % Get current position
lgdPos(2) = lgdPos(2) - 0.11; % Move legend down
lgd.Position = lgdPos; % Set new position

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
N = 60;figure,
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

%% Customized Scatter Plot with Connected Data Points Across Groups
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other')

close all
DataArray = nan(5,2);
DataArray(1:5,1) = [6,7,8,1,3];
DataArray(1:5,2) = [2,3,4,1,6];

Colors = [0.4922 0.0039 0.9063; 0.9922 0.8672 0.0039];

figure;
% Assuming 'MarkerSize' is a supported property, adjust its value to increase marker size
[xPositions, yPositions, ~, ~] = UnivarScatter(DataArray,'Label',{'meg','fmri'},'MarkerFaceColor',Colors, 'PointSize',100);
ylabel('Laterality','FontSize', 16);
xlabel('Modality','FontSize', 16);
set(gca,'color','none');
title('Word-recognition','fontsize',16)
set(gca,'FontName','HelveticaNeueLT Std Lt');
hold on

f = [xPositions, yPositions];
for j=1:length(f)
   line([f(j,1),f(j,2)],[f(j,3),f(j,4)], 'LineWidth', 1); % Adjust line width if necessary
end

% Adjust x-axis limits
xlim([0.5 2.5]);

%% Legned stop
legend('AutoUpdate', 'off')

h2 = fill(x2, inBetween, colr(j,:), 'EdgeColor', 'none', 'facealpha', 0.1);
h2.Annotation.LegendInformation.IconDisplayStyle = 'off'; % make the legend for step plot off

%% Legend outside
lgd = legend([labels]);
set(lgd, 'Box', 'off');
legend('Location', 'eastoutside');

%% Legend fit into subplots
lgd = legend(intervalTypes, 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', length(intervalTypes));
lgdPos = lgd.Position; % Get current position
lgdPos(2) = lgdPos(2) - 0.10; % Move legend down
lgd.Position = lgdPos; % Set new position

%% Adding panel labels
text(-0.2, 1.3, '(b)', 'Units', 'normalized', 'FontSize', 14);

%% replace 'banana' with 'pear'
cell_array = strrep(cell_array, 'banana', 'pear');

%% Export a figure
cfg = [];
cfg.outdir = outdir;
filename = ['meanLIs_', S_data_sel.s_tag];
cfg.filename = filename;
cfg.type = 'fig';
do_export_fig(cfg)
