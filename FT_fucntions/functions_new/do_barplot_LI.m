function do_barplot_LI(cfg)

sub_sel = cfg.sub_sel;
d_in = cfg.d_in;
L = length(d_in);

figure, bar(d_in, 0.4, 'FaceColor',[.5 .5 .5])
set(gca,'Xtick', 1:L,'XtickLabel',sub_sel);
set(gca,'FontSize',8,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   400   1000   500]);
set(gca,'color','none');
title(cfg.tit);
xlim([0,L+1])
grid