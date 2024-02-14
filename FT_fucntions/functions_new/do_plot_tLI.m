function do_plot_tLI(cfg_main)

network_sel = cfg_main.network_sel;
mLI_sub = cfg_main.mLI_sub;
wi = cfg_main.wi;
colr = cfg_main.colr;
net_sel_mutiple_label = cfg_main.net_label;
s_tag = cfg_main.S_data;
data_save_dir = cfg_main.outdir;

%%
figure,
clear LI_val
for j=1:length(network_sel)
    hold on
    plot(mLI_sub(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   500]);
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(s_tag)
ylabel('LI')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');


cfg = []; cfg.outdir = data_save_dir; 
cfg.filename = ['tLI-', s_tag];  
cfg.type = 'fig'; do_export_fig(cfg)

end