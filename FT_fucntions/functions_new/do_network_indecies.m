function [opt_idx23_L, opt_idx23_R] = do_network_indecies(cfg)

roi_network = cfg.roi_network;
idx23_L = cfg.idx23_L;
idx23_R = cfg.idx23_R;
net_sel =  cfg.net_sel;
groups_labels = cfg.groups_labels;

opt_groups_labels_num = [];
for i=1:length(net_sel)
    opt_groups_labels_num{i} = [num2str(roi_network(net_sel(i))), ': ', groups_labels{roi_network(net_sel(i))}{1}];
end
% disp(cell2table(opt_groups_labels_num'));

opt_idx23_L = []; opt_idx23_R = [];
for i=1:length(net_sel)
    opt_idx23_L = [opt_idx23_L, idx23_L{roi_network(net_sel(i))}];
    opt_idx23_R = [opt_idx23_R, idx23_R{roi_network(net_sel(i))}];
end

opt_idx23_L = sort(opt_idx23_L);
opt_idx23_R = sort(opt_idx23_R);

end