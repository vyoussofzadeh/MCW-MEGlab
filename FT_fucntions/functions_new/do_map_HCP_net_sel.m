function [idx_R_all, idx_L_all] = do_map_HCP_net_sel(cfg_main)

Data_hcp_atlas = cfg_main.Data_hcp_atlas;
idx_R = cfg_main.idx_R;
idx_L = cfg_main.idx_L;
net_sel = cfg_main.network_sel;
src_fname = cfg_main.src_fname;

%%
idx_R_all = [];
for i=1:length(net_sel)
    idx_R_all = [idx_R_all, idx_R{net_sel(i)}];
end

idx_L_all = [];
for i=1:length(net_sel)
    idx_L_all = [idx_L_all, idx_L{net_sel(i)}];
end

% idx_R_whole = []; for i=1:length(idx_R), idx_R_whole = [idx_R_whole, idx_R{i}]; end
% idx_L_whole = []; for i=1:length(idx_L), idx_L_whole = [idx_L_whole, idx_L{i}]; end

% close all
cfg = [];
cfg.atlas = Data_hcp_atlas.atlas;
cfg.src_fname = src_fname;
cfg.sel = 'whole'; % 'whole', 'left', 'right', 'roi';
cfg.index_L = idx_L_all';
cfg.index_R = idx_R_all';
cfg.rois = Data_hcp_atlas.rois;
cfg.rois_sel = 1:size(idx_L_all,2);
if isfield(cfg_main,'fixedcolor')
    cfg.fixedcolor = cfg_main.fixedcolor;
end
do_plot_HCP_atlas(cfg)

% if cfg_main.export ==1   
%     cfg = [];
%     cfg.outdir = cfg_main.savedir;
%     if
%     cfg.filename = [Data_hcp_atlas.groups_labels{cfg.network_sel}, '-hcp'];
%     cfg.type = 'fig';
%     do_export_fig(cfg)
% end

end