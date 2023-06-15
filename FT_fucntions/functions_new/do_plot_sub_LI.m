function [mm] = do_plot_sub_LI(cfg_main)

net_sel_mutiple_label = cfg_main.net_sel_mutiple_label;
LI_sub = cfg_main.LI_sub;
wi = cfg_main.wi;
subsel = cfg_main.subsel;

for i= subsel
    mm = squeeze(mean(LI_sub(:,i,:)));
    if cfg_main.plotflag == 1
        figure,
        for j=1:length(net_sel_mutiple_label)
            plot(squeeze(LI_sub(j,i,:))),
            val = round(mean(wi(:,1),2),2);
            set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
            set(gca,'FontSize',8,'XTickLabelRotation',90);
            set(gcf, 'Position', [1000   400   1100   300]);
            hold on
            title(['sub:', num2str(i)])
        end
        plot(mm,'LineWidth',3),
        legend([net_sel_mutiple_label; 'mean'])
    end
end
end