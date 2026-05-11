function plot_rows_tiled(X, fs, chs)
if nargin<3 || isempty(chs), chs = 1:min(12,size(X,1)); end
X = X(chs,:); [C,T] = size(X); t = (0:T-1)/fs;

figure('Color','w');
tiledlayout(C,1,'Padding','compact','TileSpacing','compact');
for i = 1:C
    nexttile; plot(t, X(i,:), 'k', 'LineWidth', 1);
    xlim([t(1) t(end)]); grid on;
    ylabel(sprintf('Ch %d', chs(i)));
    if i==1, title('One channel per row'); end
    if i<C, set(gca,'XTickLabel',[]); else, xlabel('Time (s)'); end
end
end
