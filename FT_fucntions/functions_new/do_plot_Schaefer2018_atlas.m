function [idx_L32, idx_R32, groups_labels_num] = do_plot_Schaefer2018_atlas(cfg)

group_members = cfg.group_members;
rois = cfg.rois;
group_labels = cfg.group_labels;
sel = cfg.roi_sel;
lat_index = cfg.lat_index;

idx_L = [];
for i=1:length(lat_index)
    grois = group_members{i};
    idx = [];
    for j=1:length(grois)
        idx(j) = strmatch(num2str(grois(j)),lat_index);
    end
    idx_L{i} = idx;
end

idx_R = [];
for i=1:length(group_members)
    grois = group_members{i}(1:end);
    idx = [];
    for j=1:length(grois)
        idx(j) = strmatch(grois{j},rois);
    end
    idx_R{i} = idx;
end

%%
idx_L = 1:200;
idx_R = 201:400;

idx_LR32 = [idx_L,idx_R];

Scouts = cfg.atlas.Scouts;
nScouts = length(Scouts);
src_fname = cfg.src_fname;
src = ft_read_headshape(src_fname);

%- Whole atlas
vertexcolor = zeros(size(src.pos,1), 3);
for iScout=1:nScouts
    index = Scouts(iScout).Vertices;
    if ~isempty(index)
        vertexcolor(index,:) = repmat(Scouts(iScout).Color,  length(index), 1);
        vertexcolor(index,:) = repmat([0.5,0.5,0.5],  length(index), 1);
    end
end


groups_labels_num = [];
for i=1:length(group_labels)
    groups_labels_num{i} = [num2str(i), ': ', group_labels{i}];
end
disp(cell2table(groups_labels_num'));
%         sel = input('enter rois:');
for iScout=1:length(sel)
    for j=1:length(idx_L{sel(iScout)})
        index = Scouts((idx_L{sel(iScout)}(j))).Vertices;
        if ~isempty(index)
            vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
        end
    end
end
disp((groups_labels_num(sel)'));

groups_labels_num = [];
for i=1:length(group_labels)
    groups_labels_num{i} = [num2str(i), ': ', group_labels{i}];
end
disp(cell2table(groups_labels_num'));
%         sel = input('enter rois:');
for iScout=1:length(sel)
    for j=1:length(idx_L{sel(iScout)})
        index = Scouts((idx_L{sel(iScout)}(j))).Vertices;
        if ~isempty(index)
            vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
        end
    end
end
disp((groups_labels_num(sel)'));


figure,
cfg = [];
cfg.view = [-180,-90;0,90;-90,0; 90,0; 0, 0];
cfg.position = [800   800   1000   300];
cfg.color = (viridis(256));
cfg.title = [''];
cfg.alpha = 1;
cfg.coor = [];
cfg.surf = src;
cfg.d_in = vertexcolor;
do_surfplot(cfg);



