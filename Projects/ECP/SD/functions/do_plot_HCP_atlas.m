function do_plot_atlas(cfg)

idx_L = cfg.lat_index(:,1);
idx_R = cfg.lat_index(:,2);
idx_lr = [idx_L;idx_R];

rois = cfg.rois;

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
        % all ROIs
        for iScout=1:length(idx_lr)
            index = Scouts(idx_lr(iScout)).Vertices;
            if ~isempty(index)
                vertexcolor(index,:) = repmat(Scouts(idx_lr(iScout)).Color,  length(index), 1);
            end
        end
    case 'left'
        
        %
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
        sel = input('enter rois (1-180):');
        idx_L_Sel = idx_L(sel);
        disp(rois(idx_L_Sel)')
        for iScout=1:length(idx_L_Sel)
            index = Scouts(idx_L_Sel(iScout)).Vertices;
            if ~isempty(index)
                vertexcolor(index,:) = repmat(Scouts(idx_L_Sel(iScout)).Color,  length(index), 1);
            end
        end
        
end
% Visualisation de l'atlas Desikan_killiany
% close all
% figure;
% ft_plot_mesh(src, 'faecolor', 'brain',  'vertexcolor', ...
%     vertexcolor, 'facealpha', 1);
% view(-3, 2);

% addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/Colormaps-from-MatPlotLib2.0')

figure
cfg = [];
cfg.view = [-180,-90;0,90;-90,0; 90,0; 0, 0];
cfg.position = [800   800   1000   300];
cfg.color = (viridis(256));
cfg.title = ['']; cfg.alpha = 1; cfg.coor = [];
cfg.surf = src;
cfg.d_in = vertexcolor;
do_surfplot(cfg);
% title([num2str(sel), ': ', rois{idx_L_Sel}])