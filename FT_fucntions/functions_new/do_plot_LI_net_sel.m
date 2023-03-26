function do_plot_LI_net_sel(cfg_main)

Data_hcp_atlas = cfg_main.Data_hcp_atlas;
S_data_sel = cfg_main.S_data_sel;
BS_data_dir = cfg_main.BS_data_dir;
idx_R = cfg_main.idx_R;
idx_L = cfg_main.idx_L;
wi = cfg_main.wi;
net_tag = cfg_main.net_tag;
network_sel = cfg_main.network_sel;

%%
cfg = [];
cfg.sinput = S_data_sel.s_avg_Input;
cfg.BS_data_dir = BS_data_dir;
cfg.atlas = Data_hcp_atlas.atlas;
cfg.thre = 0;
cfg.wi = wi;
cfg.index_L = idx_L{network_sel}; % glasser_lateral is not symmetric!
cfg.index_R = idx_R{network_sel}; % glasser_lateral is not symmetric!
cfg.fplot = 1;
cfg.tit = net_tag;
[LI, LI_max] = do_lat_analysis_asymetric(cfg);

end