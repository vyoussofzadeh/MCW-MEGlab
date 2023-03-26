function do_plot_pow_net_all(cfg_main)

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
lgd = []; k=1;
for i=1:length(Data_hcp_atlas.groups_labels)
    lgd{k} = [Data_hcp_atlas.groups_labels{i}, '-R'];
    k=1+k;
    lgd{k} = [Data_hcp_atlas.groups_labels{i}, '-L'];
    k=1+k;
end

%%
figure,
for i=1:length(idx_L)
    cfg = []; cfg.sinput = S_data_sel.s_avg_Input;
    cfg.BS_data_dir = BS_data_dir;
    cfg.atlas = Data_hcp_atlas.atlas; cfg.thre = thre; cfg.wi = wi;
    cfg.index_L = idx_L{i}; % glasser_lateral is not symmetric!
    cfg.index_R = idx_R{i}; % glasser_lateral is not symmetric!
    cfg.fplot = 0;
    cfg.tit = [''];
    [pow_left, pow_right] = do_pow_analysis_asymetric(cfg);
    pow_left_norm = normalize(pow_left, 2);
    pow_right_norm = normalize(pow_right, 2);
    hold on
    plot(pow_right_norm, 'Color', colr(i,:)),
    plot(pow_left_norm, 'Color', colr(i,:), 'LineWidth', 2),
end
set(gcf, 'Position', [1000   800   1000   700]);
title(S_data_sel.s_tag),
val = round(mean(wi(:,1),2),2);
lgnd = legend(lgd);
set(lgnd,'color','none');
set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
set(gca,'FontSize',8,'XTickLabelRotation',90);
xlabel('temporal windows (sec)')
ylabel('Power')
set(gca,'color','none');

end