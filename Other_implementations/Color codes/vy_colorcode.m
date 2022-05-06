colr = hsv(nScouts);
colr = hot(nScouts);
% colr = viridis(nScouts);
addpath('/data/MEG/Vahab/Github/MCW-MEGlab/tools/helpful_tools');
colr = viridis(nScouts);

% colr = hsv(nScouts);
addpath('/data/MEG/Vahab/Github/MCW-MEGlab/FT/functions/External/brewermap');
colr = distinguishable_colors(nScouts);

%% Example
colr = distinguishable_colors(9);
val = ...;
figure, 
set(gca, 'ColorOrder', colr, 'NextPlot', 'replacechildren');
plot(val),
legend('blahblah..')