function do_plot_LI_net_all(cfg_main)

Data_hcp_atlas = cfg_main.Data_hcp_atlas;
S_data_sel = cfg_main.S_data_sel;
BS_data_dir = cfg_main.BS_data_dir;
idx_R = cfg_main.idx_R;
idx_L = cfg_main.idx_L;
wi = cfg_main.wi;
% atlas = cfg_main.atlas;
% net_tag = cfg_main.net_tag;
thre = cfg_main.thre;


%%
colr = distinguishable_colors(length(idx_L));

%%
figure,
for i=1:length(idx_L)
    cfg = []; cfg.sinput = S_data_sel.s_avg_Input; 
    cfg.BS_data_dir = BS_data_dir;
    cfg.atlas = Data_hcp_atlas.atlas; 
    cfg.thre = thre; 
    cfg.wi = wi;
    cfg.index_L = idx_L{i}; % glasser_lateral is not symmetric!
    cfg.index_R = idx_R{i}; % glasser_lateral is not symmetric!
    cfg.fplot = 0; 
    cfg.tit = [''];
    [LI, ~] = do_lat_analysis_asymetric(cfg);
    hold on
    plot(LI, 'Color', colr(i,:)),
end
set(gcf, 'Position', [1000   400   1000   300]);
title(S_data_sel.s_tag),
val = round(mean(wi(:,1),2),2);
lgnd = legend(Data_hcp_atlas.groups_labels);
set(lgnd,'color','none');
set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
set(gca,'FontSize',8,'XTickLabelRotation',90);
xlabel('temporal windows (sec)')
ylabel('LI')
set(gca,'color','none');

end