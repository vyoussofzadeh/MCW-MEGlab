function vy_plotrbar(cfg, datain)


bar(datain, 0.4), title(cfg.titletag)
set(gca,'Xtick', 1:2,'XtickLabel',cfg.XtickLabeltag);
xlabel(cfg.xlabeltag), ylabel(cfg.ylabeltag),
% set(gca,'FontSize',10,'XTickLabelRotation',90);
% grid on
set(gca,'color','none');

