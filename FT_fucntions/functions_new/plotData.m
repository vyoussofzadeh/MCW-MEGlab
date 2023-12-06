function plotData(cfg)

data = cfg.data;
labels = cfg.labels;
colors = cfg.colors;
titleText = cfg.titleText;
wi = cfg.wi;

figure,
for j=1:length(labels)
    hold on
    plot(mean(wi'),data(j, :), 'LineWidth', 3, 'color', colors(j, :));
end
set(gca,'color','none');
legend(labels);

if length(labels) > 5
legend(labels,'Location','southoutside', 'NumColumns', 5)
end

title(titleText);
ylabel('LI');
xlabel('time');
set(gca, 'color', 'none');

end
