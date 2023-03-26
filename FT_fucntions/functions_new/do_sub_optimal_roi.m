function wi_sub_max = do_sub_optimal_roi(cfg_main)

Data_hcp_atlas = cfg_main.Data_hcp_atlas;
S_data_sel = cfg_main.S_data_sel;
wi_sub_max = cfg_main.wi_sub_max;
thre = cfg_main.thre;
BS_data_dir = cfg_main.BS_data_dir;
idx_L_whole = cfg_main.idx_L_whole;
idx_R_whole = cfg_main.idx_R_whole;

cfg = [];
cfg.sinput = S_data_sel.sFiles_in;
cfg.BS_data_dir = BS_data_dir;
cfg.atlas = Data_hcp_atlas.atlas; cfg.thre = 0; cfg.fplot = 0;
% cfg.lat_index = [opt_idx_L; opt_idx_R]';
% cfg.index_L = idx_L_all;
% cfg.index_R = idx_R_all;
cfg.index_L = idx_L_whole;
cfg.index_R = idx_R_whole;


sFiles_in = S_data_sel.sFiles_in;

ft_progress('init', 'text',     'please wait ...');
clear m_LI_max_sub LI_sub
for i=1:length(sFiles_in)
    ft_progress(i/length(sFiles_in), 'Processing subjects %d from %d', i, length(sFiles_in));
    pause(0.1);
    cfg.sinput = sFiles_in{i};
    cfg.wi = wi_sub_max(i,:);
    [LI, wi_max] = do_lat_analysis_asymetric(cfg);
    LI_sub{i} = LI;
    m_LI_max_sub(i) = mean(LI);
    %     wi_sub_max(i,:) = wi_max;
end
ft_progress('close')

% Plot LIs
cfg = [];
cfg.sub_sel = S_data_sel.sFiles_subid;
cfg.d_in = m_LI_max_sub;
cfg.tit = ['LIs (network): ', S_data_sel.taskcon, '-', S_data_sel.subcon];
do_barplot_LI(cfg)
set(gcf, 'Position', [1000   400   1000   300]);

% figure, bar(m_LI_sub, 0.4)
% set(gca,'Xtick', 1:length(sFiles_in),'XtickLabel',sub);
% set(gca,'FontSize',8,'XTickLabelRotation',90);
% set(gcf, 'Position', [1000   400   1000   500]);
% set(gca,'color','none');
% disp(wi_sub_max)

mean(wi_sub_max)
std(wi_sub_max)

end