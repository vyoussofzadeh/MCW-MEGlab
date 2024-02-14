function [mLI_sub_pt] = do_tLI_diff(cfg_main)

d1 = cfg_main.d1;
d2 = cfg_main.d2;
network_sel = cfg_main.spec.network_sel;
net_sel_mutiple_label = cfg_main.spec.net_label;
s_tag = cfg_main.spec.s_tag;
wi = cfg_main.spec.wi;
data_save_dir = cfg_main.spec.outdir;

%%

mLI_sub1 = squeeze(mean(d1,2)); 
mLI_sub1 = [mLI_sub1; mean(mLI_sub1([1,2,6],:),1)];

mLI_sub2 = squeeze(mean(d2,2)); 
mLI_sub2 = [mLI_sub2; mean(mLI_sub2([1,2,6],:),1)];

mLI_sub_pt = mLI_sub1 - mLI_sub2;

%%

cfg = []; 
cfg.network_sel = network_sel; %[network_sel, 9]; 
cfg.wi = wi; 
cfg.colr = distinguishable_colors(length(cfg.network_sel)); 
cfg.net_label = [net_sel_mutiple_label; 'Ang-Fron-Temp']; 
cfg.outdir = data_save_dir;
cfg.mLI_sub = mLI_sub_pt; 
cfg.S_data = s_tag; 
do_plot_tLI(cfg)