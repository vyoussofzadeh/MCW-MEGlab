function [idx_R_whole, idx_L_whole]  = do_LI_avg(cfg_main)

idx_R = cfg_main.idx_R;
idx_L = cfg_main.idx_L;
S_data_sel = cfg_main.S_data_sel;
BS_data_dir = cfg_main.BS_data_dir;
Data_hcp_atlas = cfg_main.Data_hcp_atlas;
wi = cfg_main.wi;


%%
idx_R_whole = []; for i=1:length(idx_R), idx_R_whole = [idx_R_whole, idx_R{i}]; end
idx_L_whole = []; for i=1:length(idx_L), idx_L_whole = [idx_L_whole, idx_L{i}]; end

cfg = [];
cfg.sinput = S_data_sel.s_avg_Input;
cfg.BS_data_dir = BS_data_dir;
cfg.atlas = Data_hcp_atlas.atlas;
cfg.thre = cfg_main.thre;
cfg.wi = wi;
cfg.index_L = idx_L_whole; % glasser_lateral is not symmetric!
cfg.index_R = idx_R_whole; % glasser_lateral is not symmetric!
cfg.fplot = 1;
cfg.tit = 'Whole-brain, avg';%['Whole-brain, avg:', s_tag];
do_lat_analysis_asymetric(cfg);

end