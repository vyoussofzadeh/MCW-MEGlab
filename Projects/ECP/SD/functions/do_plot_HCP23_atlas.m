function do_plot_HCP23_atlas(cfg)

groups = cfg.groups;
rois = cfg.rois;
groups_labels = cfg.groups_labels;
% idx_lr = cfg.idx_lr;

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
% length(idx_lr32)


% idx_L = cfg.lat_index(:,1);
% idx_R = cfg.lat_index(:,2);
% idx_lr = [idx_L;idx_R];


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
        %         for iScout=1:length(idx_lr)
        %             index = Scouts(idx_lr(iScout)).Vertices;
        %             if ~isempty(index)
        %                 vertexcolor(index,:) = repmat(Scouts(idx_lr(iScout)).Color,  length(index), 1);
        %             end
        %         end
        
        for iScout=1:length(idx_LR32)
            for j=1:length(idx_LR32{iScout})
                index = Scouts((idx_LR32{iScout}(j))).Vertices;
                if ~isempty(index)
                    vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
                end
            end
        end
        
    case 'left'
        
        %
        % left ROIs
        for iScout=1:length(idx_L32)
            for j=1:length(idx_L32{iScout})
                index = Scouts((idx_L32{iScout}(j))).Vertices;
                if ~isempty(index)
                    vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
                end
            end
        end
%         for iScout=1:length(idx_L)
%             index = Scouts(idx_L(iScout)).Vertices;
%             if ~isempty(index)
%                 vertexcolor(index,:) = repmat(Scouts(idx_L(iScout)).Color,  length(index), 1);
%             end
%         end
        
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
%         for iScout=1:length(idx_R)
%             index = Scouts(idx_R(iScout)).Vertices;
%             if ~isempty(index)
%                 vertexcolor(index,:) = repmat(Scouts(idx_R(iScout)).Color,  length(index), 1);
%             end
%         end
        
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
%                     vertexcolor(index,:) = repmat(Scouts(idx_lr(iScout)).Color,  length(index), 1);
                    vertexcolor(index,:) = repmat(Scouts((iScout)).Color,  length(index), 1);
                end
            end
        end
        
        
        
%         sel = input('enter rois (1-180):');
%         idx_L_Sel = idx_L(sel);
%         disp(rois(idx_L_Sel)')
%         for iScout=1:length(idx_L_Sel)
%             index = Scouts(idx_L_Sel(iScout)).Vertices;
%             if ~isempty(index)
%                 vertexcolor(index,:) = repmat(Scouts(idx_L_Sel(iScout)).Color,  length(index), 1);
%             end
%         end
        
end
% Visualisation de l'atlas Desikan_killiany
% close all
% figure;
% ft_plot_mesh(src, 'faecolor', 'brain',  'vertexcolor', ...
%     vertexcolor, 'facealpha', 1);
% view(-3, 2);

% addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/Colormaps-from-MatPlotLib2.0')

% figure
% cfg = [];
% cfg.view = [-180,-90;0,90;-90,0; 90,0; 0, 0];
% cfg.position = [800   800   1000   300];
% cfg.color = (viridis(256));
% cfg.title = ['']; cfg.alpha = 1; cfg.coor = [];
% do_surfplot(cfg,src,vertexcolor);
% % title([num2str(sel), ': ', rois{idx_L_Sel}])

figure
cfg = [];
cfg.view = [-180,-90;0,90;-90,0; 90,0; 0, 0];
cfg.position = [800   800   1000   300];
cfg.color = (viridis(256));
cfg.title = ['']; cfg.alpha = 1; cfg.coor = [];
do_surfplot(cfg,src,vertexcolor);


%%
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
% 
% %%
% idx_LR32 = [idx_L32,idx_R32];
% length(idx_lr32)
% 
% 
% index = [];
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
% 
% % left ROIs
% % for iScout=1:length(idx_L32)
% %     for j=1:length(idx_L32{iScout})
% %         index = Scouts((idx_L32{iScout}(j))).Vertices;
% %         if ~isempty(index)
% %             vertexcolor(index,:) = repmat(Scouts(idx_lr(iScout)).Color,  length(index), 1);
% %         end
% %     end
% % end
% 
% % % right ROIs
% % for iScout=1:length(idx_R32)
% %     for j=1:length(idx_R32{iScout})
% %         index = Scouts((idx_R32{iScout}(j))).Vertices;
% %         if ~isempty(index)
% %             vertexcolor(index,:) = repmat(Scouts(idx_lr(iScout)).Color,  length(index), 1);
% %         end
% %     end
% % end
% 
% 
% % close all
% % left_sel ROIs
% groups_labels_num = [];
% for i=1:length(groups_labels)
%     groups_labels_num{i} = [num2str(i), ': ', groups_labels{i}{1}];
% end
% disp(cell2table(groups_labels_num'));
% sel = input('enter rois:');
% for iScout=1:length(sel)
%     for j=1:length(idx_L32{sel(iScout)})
%         index = Scouts((idx_L32{sel(iScout)}(j))).Vertices;
%         if ~isempty(index)
%             vertexcolor(index,:) = repmat(Scouts(idx_lr(iScout)).Color,  length(index), 1);
%         end
%     end
% end
% 
% 
% % Visualisation de l'atlas Desikan_killiany
% % close all
% % figure;
% % ft_plot_mesh(src, 'faecolor', 'brain',  'vertexcolor', ...
% %     vertexcolor, 'facealpha', 1);
% % view(-3, 2);
% 
% addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/Colormaps-from-MatPlotLib2.0')
% 
% figure
% cfg = [];
% cfg.view = [-180,-90;0,90;-90,0; 90,0; 0, 0];
% cfg.position = [800   800   1000   300];
% cfg.color = (viridis(256));
% cfg.title = ['']; cfg.alpha = 1; cfg.coor = [];
% do_surfplot(cfg,src,vertexcolor);
% % title([num2str(sel), ': ', rois{idx_L_Sel}])