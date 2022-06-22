% pause, close all,
cfg           = [];
cfg.method    = 'degrees';
cfg.parameter = par;
cfg.threshold = .5;
net = ft_networkanalysis(cfg,source_conn1);

f_range = round(net.freq);
f_idx = 1:2:length(net.freq);
figure, imagesc(net.degrees),
set(gca,'Xtick', f_idx,'XtickLabel',f_range(f_idx));
set(gca,'Ytick', 1:size(net.degrees,1),'YtickLabel',1:size(net.degrees,1));
set(gca,'FontSize',9,'XTickLabelRotation',90);
set(gcf, 'Position', [800   400   500   1300]);

for i=1:length(vs_roi1.label)
    roi_label{i} = [num2str(i), ':', vs_roi1.label{i}];
end
disp(roi_label');

figure, plot(mean(net.degrees)),
set(gca,'Xtick', f_idx,'XtickLabel',f_range(f_idx));
set(gca,'Ytick', 1:size(net.degrees,1),'YtickLabel',1:size(net.degrees,1));
set(gca,'FontSize',9,'XTickLabelRotation',90);