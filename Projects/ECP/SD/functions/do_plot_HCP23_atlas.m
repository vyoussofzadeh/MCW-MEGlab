function [idx_L32, idx_R32] = do_plot_HCP23_atlas(cfg)

groups = cfg.groups;
rois = cfg.rois;
groups_labels = cfg.groups_labels;

idx_L32 = [];
for i=1:length(groups)
    grois = groups{i}(2:end);
    idx = [];
    for j=1:length(grois)
        idx(j) = strmatch(['L_',grois{j}, '_ROI'],rois);
    end
    idx_L32{i} = idx;
end

idx_R32 = [];
for i=1:length(groups)
    grois = groups{i}(2:end);
    idx = [];
    for j=1:length(grois)
        idx(j) = strmatch(['R_',grois{j}, '_ROI'],rois);
    end
    idx_R32{i} = idx;
end

%%
idx_LR32 = [idx_L32,idx_R32];

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

switch cfg.sel
    
    case 'whole'
        for iScout=1:length(idx_LR32)
            for j=1:length(idx_LR32{iScout})
                index = Scouts((idx_LR32{iScout}(j))).Vertices;
                if ~isempty(index)
                    vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
                end
            end
        end
        
    case 'left'
        
        % left ROIs
        for iScout=1:length(idx_L32)
            for j=1:length(idx_L32{iScout})
                index = Scouts((idx_L32{iScout}(j))).Vertices;
                if ~isempty(index)
                    vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
                end
            end
        end
        
    case 'right'
        % % right ROIs
        for iScout=1:length(idx_R32)
            for j=1:length(idx_R32{iScout})
                index = Scouts((idx_R32{iScout}(j))).Vertices;
                if ~isempty(index)
                    vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
                end
            end
        end
        
    case 'roi'
        
        groups_labels_num = [];
        for i=1:length(groups_labels)
            groups_labels_num{i} = [num2str(i), ': ', groups_labels{i}{1}];
        end
        disp(cell2table(groups_labels_num'));
        sel = input('enter rois:');
        for iScout=1:length(sel)
            for j=1:length(idx_L32{sel(iScout)})
                index = Scouts((idx_L32{sel(iScout)}(j))).Vertices;
                if ~isempty(index)
                    vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
                end
            end
        end
        disp((groups_labels_num(sel)'));
end


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



