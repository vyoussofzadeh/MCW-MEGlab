function [opt_idx_L, opt_idx_R] = do_network_indecies2(cfg)

roi_network = cfg.roi_network;
idx_L = cfg.idx_L;
idx_R = cfg.idx_R;
net_sel =  cfg.net_sel;
groups_labels = cfg.groups_labels;

opt_groups_labels_num = [];
for i=1:length(net_sel)
    opt_groups_labels_num{i} = [num2str(roi_network(net_sel(i))), ': ', groups_labels{roi_network(net_sel(i))}];
end
% disp(cell2table(opt_groups_labels_num'));

opt_idx_L = []; opt_idx_R = [];
for i=1:length(net_sel)
    opt_idx_L = [opt_idx_L, idx_L{roi_network(net_sel(i))}];
    opt_idx_R = [opt_idx_R, idx_R{roi_network(net_sel(i))}];
end

opt_idx_L = sort(opt_idx_L);
opt_idx_R = sort(opt_idx_R);

end