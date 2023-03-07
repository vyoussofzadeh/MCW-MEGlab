function do_plot_HCP_atlas(cfg_main)

idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;

idx_lr = [idx_L;idx_R];

if size(idx_lr,2) > 1
    error('check the size of index')
end

rois = cfg_main.rois;

Scouts = cfg_main.atlas.Scouts;
nScouts = length(Scouts);
src_fname = cfg_main.src_fname;
src = ft_read_headshape(src_fname);

% addpath('/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master/template/anatomy')
% surface_pial_both = load('surface_inflated_both.mat');
% src1 = surface_pial_both.mesh;

%- Whole atlas
vertexcolor = zeros(size(src.pos,1), 3);
for iScout=1:nScouts
    index = Scouts(iScout).Vertices;
    if ~isempty(index)
        vertexcolor(index,:) = repmat(Scouts(iScout).Color,  length(index), 1);
        vertexcolor(index,:) = repmat([0.5,0.5,0.5],  length(index), 1);
    end
end

switch cfg_main.sel
    
    case 'whole'
        % all ROIs
        for iScout=1:length(idx_lr)
            index = Scouts(idx_lr(iScout)).Vertices;
            if ~isempty(index)
                vertexcolor(index,:) = repmat(Scouts(idx_lr(iScout)).Color,  length(index), 1);
            end
        end
    case 'left'
        % left ROIs
        for iScout=1:length(idx_L)
            index = Scouts(idx_L(iScout)).Vertices;
            if ~isempty(index)
                vertexcolor(index,:) = repmat(Scouts(idx_L(iScout)).Color,  length(index), 1);
            end
        end
        
    case 'right'
        % % right ROIs
        for iScout=1:length(idx_R)
            index = Scouts(idx_R(iScout)).Vertices;
            if ~isempty(index)
                vertexcolor(index,:) = repmat(Scouts(idx_R(iScout)).Color,  length(index), 1);
            end
        end
        
    case 'roi'
        % % close all
        % % left_sel ROIs
        sel = cfg_main.rois_sel; %input('enter rois (1-180):');
        idx_L_Sel = idx_L(sel);
        disp(rois(idx_L_Sel)')
        
        %         idx_R_Sel = idx_R(sel);
        %         disp(rois(idx_R_Sel)')
        
        
        for iScout=1:length(idx_L_Sel)
            index = Scouts(idx_L_Sel(iScout)).Vertices;
            if ~isempty(index)
                vertexcolor(index,:) = repmat(Scouts(idx_L_Sel(iScout)).Color,  length(index), 1);
            end
        end
        
        %         for iScout=1:length(idx_R_Sel)
        %             index = Scouts(idx_R_Sel(iScout)).Vertices;
        %             if ~isempty(index)
        %                 vertexcolor(index,:) = repmat(Scouts(idx_R_Sel(iScout)).Color,  length(index), 1);
        %             end
        %         end
        
end

switch cfg_main.sel
    case 'roi'
        %- Left
        src_L = src;
        src_L.tri = src_L.tri(1:14980,:);
        
        figure
        cfg = [];
        cfg.view = [-180,-90; 0,90;-90,0; 90,0;];
        cfg.position = [800   800   900   200];
        cfg.color = (viridis(256));
        cfg.title = ['LH: roi', cfg_main.title];
        cfg.alpha = 1; cfg.coor = [];
        cfg.surf = src_L;
        cfg.d_in = vertexcolor;
        do_surfplot(cfg);
        
        %- Right
        src_R = src;
        src_R.tri = src_R.tri(14981:end,:);
        
        figure
        cfg = [];
        cfg.view = [-180,-90; 0,90;-90,0; 90,0;];
        cfg.position = [800   500   900   200];
        cfg.color = (viridis(256));
        cfg.title = ['RH: roi', cfg_main.title];
        cfg.alpha = 1; cfg.coor = [];
        cfg.surf = src_R;
        cfg.d_in = vertexcolor;
        do_surfplot(cfg);
        
    otherwise
        figure
        cfg = [];
        cfg.view = [-180,-90;0,90;-90,0; 90,0; 0, 0];
        cfg.position = [800   800   1000   300];
        cfg.color = (viridis(256));
        cfg.title = [''];
        cfg.alpha = 1; cfg.coor = [];
        cfg.surf = src;
        cfg.d_in = vertexcolor;
        do_surfplot(cfg);        
end


