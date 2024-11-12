function [idx_L, idx_R, src]  = do_plot_hcp_network2(cfg_main)

clc
% src_fname = cfg_main.src_fname;
Data_hcp_atlas = cfg_main.Data_hcp_atlas;
network_sel = cfg_main.network_sel;

% cfg = [];
% cfg.atlas = Data_hcp_atlas.atlas;
% cfg.src_fname = src_fname;
% cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
% cfg.rois = Data_hcp_atlas.rois;
% cfg.group_labels = Data_hcp_atlas.groups_labels;
% cfg.group_members = Data_hcp_atlas.glass_net_L_label;
% cfg.roi_sel = network_sel;
% do_plot_HCP6_atlas(cfg);
% cfg.group_members = Data_hcp_atlas.glass_net_R_label;
% do_plot_HCP6_atlas(cfg);
% cfg.sel = 'whole'; % 'whole', 'left', 'right', 'roi';
% cfg.group_members = [glass_net_R_label,glass_net_L_label];
% do_plot_HCP6_atlas(cfg);

% close all
cfg = [];
cfg.atlas = Data_hcp_atlas.atlas;
% cfg.src_fname = src_fname;
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.rois = Data_hcp_atlas.rois;
cfg.group_labels = Data_hcp_atlas.groups_labels;
cfg.group_members = Data_hcp_atlas.glass_net_L_label;
cfg.roi_sel = network_sel;
cfg.plotflag = cfg_main.plotflag;
if isfield(cfg_main,'fixedcolor')
    cfg.fixedcolor = cfg_main.fixedcolor;
end
% net_sel = [6];
% cfg.roi_sel = net_sel;
cfg.group_members = Data_hcp_atlas.glass_net_L_label;
[idx_L, ~, ~, src] = do_plot_HCP7_atlas(cfg);
cfg.group_members = Data_hcp_atlas.glass_net_R_label;
[idx_R, ~, ~] = do_plot_HCP7_atlas(cfg);

end