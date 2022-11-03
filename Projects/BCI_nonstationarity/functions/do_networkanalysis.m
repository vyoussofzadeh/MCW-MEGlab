function net = do_networkanalysis(mcfg, source_conn1)

cfg           = [];
cfg.method    = 'degrees';
cfg.parameter = mcfg.par;
cfg.threshold = .5;
net = ft_networkanalysis(cfg,source_conn1);

f_range = round(net.freq);
f_idx = 1:2:length(net.freq);
figure, imagesc(net.degrees),
set(gca,'Xtick', f_idx,'XtickLabel',f_range(f_idx));
set(gca,'Ytick', 1:size(net.degrees,1),'YtickLabel',1:size(net.degrees,1));
set(gca,'FontSize',9,'XTickLabelRotation',90);
set(gcf, 'Position', [800   400   500   1300]);

for i=1:length(mcfg.label)
    roi_label{i} = [num2str(i), ':', mcfg.label{i}];
end
disp(roi_label');

figure, plot(mean(net.degrees)),
set(gca,'Xtick', f_idx,'XtickLabel',f_range(f_idx));
set(gca,'Ytick', 1:size(net.degrees,1),'YtickLabel',1:size(net.degrees,1));
set(gca,'FontSize',9,'XTickLabelRotation',90);