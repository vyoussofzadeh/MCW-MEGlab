function [LI_sub, m_LI_sub, wi_sub_max] = do_sub_LI(cfg_main)

Data_hcp_atlas = cfg_main.Data_hcp_atlas;
S_data_sel = cfg_main.S_data_sel;
wi = cfg_main.wi;
thre = cfg_main.thre;
BS_data_dir = cfg_main.BS_data_dir;
idx_L_whole = cfg_main.idx_L_whole;
idx_R_whole = cfg_main.idx_R_whole;

%%
cfg = [];
cfg.sinput = S_data_sel.sFiles_in;
cfg.BS_data_dir = BS_data_dir;
cfg.atlas = Data_hcp_atlas.atlas;
cfg.thre = thre;
cfg.fplot = 0;
cfg.wi = wi;
cfg.index_L = idx_L_whole;
cfg.index_R = idx_R_whole;

ft_progress('init', 'text',     'please wait ...');

sFiles_in = S_data_sel.sFiles_in;
clear m_LI_sub LI_sub wi_sub_max
for i=1:length(sFiles_in)
    ft_progress(i/length(sFiles_in), 'Processing subjects %d from %d', i, length(sFiles_in));
    pause(0.1);
    cfg.sinput = sFiles_in{i};
    cfg.tit = ['LI, sub:', num2str(i)];
    [LI, wi_max] = do_lat_analysis_asymetric(cfg);
    LI_sub{i} = LI;
    m_LI_sub(i) = nanmean(LI);
    wi_sub_max(i,:) = wi_max;
end
ft_progress('close')

% Plot LIs
cfg = [];
cfg.sub_sel = S_data_sel.sFiles_subid;
cfg.d_in = m_LI_sub;
cfg.tit = ['LIs (network): ', S_data_sel.taskcon, '-', S_data_sel.subcon];
do_barplot_LI(cfg)
set(gcf, 'Position', [1000   400   1000   300]);


figure,
for i=1:length(LI_sub)
    plot(LI_sub{i}),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    hold on
end
plot(squeeze(mean(cell2mat(LI_sub'),1)),'LineWidth',3),
set(gcf, 'Position', [1000   400   1000   300]);

end