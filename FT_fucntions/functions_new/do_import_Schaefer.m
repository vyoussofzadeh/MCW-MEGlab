function [idx_L32, idx_R32, groups_labels_num] = do_import_Schaefer()

atlas_bs_mat  = ('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/Schaefer2018/schaefer/scout_Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm_400.mat');
schaefer_bs = load(atlas_bs_mat);

% atlas_nii = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/Schaefer2018/schaefer/MNI/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii.gz';
% atlas = ft_read_atlas(atlas_nii)

for i=1:length(schaefer_bs.Scouts)
    rois{i} = schaefer_bs.Scouts(i).Label;
    region{i} = schaefer_bs.Scouts(i).Region;
    region_lr{i} = schaefer_bs.Scouts(i).Region(1);
end

%%
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = "Networks_LH_Vis_1";
opts.VariableTypes = "string";
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
Schaefer2018400Parcels7Networksorderinfo = readtable("/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/Schaefer2018/schaefer/HCP/fslr32k/gifti/Schaefer2018_400Parcels_7Networks_order_info.txt", opts);

%%
net = Schaefer2018400Parcels7Networksorderinfo.Networks_LH_Vis_1;

net_id = []; Lat_id = [];
k = 1;
for i = 1:2:length(net)  
    tmp = net{i}; 
    tkz = tokenize(tmp,'_');
    net_id{k} = tkz{3}; 
    Lat_id{k} = tkz{2}; k = k+1;
end

snet_id = cell2table(net_id');

net_list = unique(net_id);

for i=1:length(net_list)
    net_list_idx{i} = find(contains(snet_id.Var1,net_list{i})==1);
end

%%
% idx_L = 1:200;
% idx_R = 201:400;

cfg = [];
cfg.atlas = schaefer_bs;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.lat_index = Lat_id;
cfg.rois = rois;
% cfg.groups_labels = net_list;
% cfg.groups = net_list_idx;
cfg.group_labels = net_list;
cfg.group_members = net_list_idx;
% cfg.roi_sel = [12,21]; do_plot_HCP23_atlas(cfg);
% cfg.roi_sel = [13,14,15]; do_plot_HCP23_atlas(cfg);
% cfg.roi_sel = [16,17]; do_plot_HCP23_atlas(cfg);
% cfg.roi_sel = [12,21,13,14,15]; do_plot_HCP23_atlas(cfg);
% cfg.roi_sel = [11,12,13,14,16,17,20,21]; do_plot_HCP23_atlas(cfg);
cfg.roi_sel = [1:2]; do_plot_Schaefer2018_atlas(cfg);


%%

% groups = cfg.groups;
% rois = cfg.rois;
% groups_labels = cfg.groups_labels;
% sel = cfg.roi_sel;
% 
% idx_L32 = [];
% for i=1:length(groups)
%     grois = groups{i}(2:end);
%     idx = [];
%     for j=1:length(grois)
%         idx(j) = strmatch(['L_',grois{j}, '_ROI'],rois);
%     end
%     idx_L32{i} = idx;
% end
% 
% idx_R32 = [];
% for i=1:length(groups)
%     grois = groups{i}(2:end);
%     idx = [];
%     for j=1:length(grois)
%         idx(j) = strmatch(['R_',grois{j}, '_ROI'],rois);
%     end
%     idx_R32{i} = idx;
% end
% 
% %%
% idx_LR32 = [idx_L32,idx_R32];
% 
% Scouts = cfg.atlas.Scouts;
% nScouts = length(Scouts);
% src_fname = cfg.src_fname;
% src = ft_read_headshape(src_fname);
% 
% %- Whole atlas
% vertexcolor = zeros(size(src.pos,1), 3);
% for iScout=1:nScouts
%     index = Scouts(iScout).Vertices;
%     if ~isempty(index)
%         vertexcolor(index,:) = repmat(Scouts(iScout).Color,  length(index), 1);
%         vertexcolor(index,:) = repmat([0.5,0.5,0.5],  length(index), 1);
%     end
% end
% 
% switch cfg.sel
%     
%     case 'whole'
%         for iScout=1:length(idx_LR32)
%             for j=1:length(idx_LR32{iScout})
%                 index = Scouts((idx_LR32{iScout}(j))).Vertices;
%                 if ~isempty(index)
%                     vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
%                 end
%             end
%         end
%         
%     case 'left'
%         
%         % left ROIs
%         for iScout=1:length(idx_L32)
%             for j=1:length(idx_L32{iScout})
%                 index = Scouts((idx_L32{iScout}(j))).Vertices;
%                 if ~isempty(index)
%                     vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
%                 end
%             end
%         end
%         
%     case 'right'
%         % % right ROIs
%         for iScout=1:length(idx_R32)
%             for j=1:length(idx_R32{iScout})
%                 index = Scouts((idx_R32{iScout}(j))).Vertices;
%                 if ~isempty(index)
%                     vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
%                 end
%             end
%         end
%         
%     case 'roi'
%         
%         groups_labels_num = [];
%         for i=1:length(groups_labels)
%             groups_labels_num{i} = [num2str(i), ': ', groups_labels{i}{1}];
%         end
%         disp(cell2table(groups_labels_num'));
% %         sel = input('enter rois:');
%         for iScout=1:length(sel)
%             for j=1:length(idx_L32{sel(iScout)})
%                 index = Scouts((idx_L32{sel(iScout)}(j))).Vertices;
%                 if ~isempty(index)
%                     vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
%                 end
%             end
%         end
%         disp((groups_labels_num(sel)'));
% end
% 
% 
% figure,
% cfg = [];
% cfg.view = [-180,-90;0,90;-90,0; 90,0; 0, 0];
% cfg.position = [800   800   1000   300];
% cfg.color = (viridis(256));
% cfg.title = ['']; 
% cfg.alpha = 1; 
% cfg.coor = [];
% cfg.surf = src; 
% cfg.d_in = vertexcolor;
% do_surfplot(cfg);



