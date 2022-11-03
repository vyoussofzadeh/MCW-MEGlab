function [vs, vs_roi1] = do_extractvirtualsensor(mcfg, source_active)

% source_whole_temp = source_whole;
% source_whole_temp.pos     = template_grid.pos;
% source_whole_temp.dim     = template_grid.dim;
% source_whole_temp.inside  = template_grid.inside;

source_active_temp = source_active;
% source_active_temp.pos     = template_grid.pos;
% source_active_temp.dim     = template_grid.dim;
% source_active_temp.inside  = template_grid.inside;

sourcemodel = mcfg.individual_grid;

%%
% template_grid.coordsys = 'mni';
atlas = mcfg.atlas;
cfg = [];
cfg.atlas      = atlas;
cfg.roi        = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
cfg.inputcoord = 'mni';
mask           = ft_volumelookup(cfg, sourcemodel);

%%
Label = zeros(length(sourcemodel.pos),1);

for jj=1:length(atlas.tissuelabel)
    
    cfg = [];
    cfg.atlas      = atlas;
    cfg.roi        = atlas.tissuelabel{jj};  % here you can also specify a single label, i.e. single ROI
    cfg.inputcoord = 'mni';
    mask1           = ft_volumelookup(cfg, sourcemodel);
    
    
    sourcemodel3 = sourcemodel;
    sourcemodel3.inside = false(sourcemodel3.dim);
    sourcemodel3.inside(mask1==1) = true;
    
    grid2 = sourcemodel3;
    grid2.inside = reshape(grid2.inside,length(grid2.pos),1);
    
    find(grid2.inside==1)
    Label(grid2.inside==1,:)= jj;
    
end
% figure, bar(Label(Label~=0))

%% Masking using aal atlas
% pause, close all,
sourcemodel2 = sourcemodel;
sourcemodel2.inside = false(sourcemodel2.dim);
sourcemodel2.inside(mask==1) = true;

% figure; hold on;
% ft_plot_vol(vol, 'edgecolor', 'none', 'facealpha', 0.4);
% ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
grid = sourcemodel2;
grid.inside = reshape(grid.inside,length(grid.pos),1);

%% VS extraction template positions
% pause, close all,

active_nodes = find(source_active_temp.inside==1);
% create the trial field from source anal output
vs = [];
node = [];
for i = 1:length(source_active_temp.trial) % i = number trials
    for x = 1:length(active_nodes) % x = number nodes
        node(x,:) = source_active_temp.trial(i).mom{active_nodes(x)};
    end
    vs.trial{i} = node;
end

% set up labels
label_temp = num2str(active_nodes);
label_temp = cellstr(label_temp);
vs.label = label_temp;

% set up time component
for i = 1:length(source_active_temp.trial)
    vs.time{1, i} = source_active_temp.time;
end

vs.pos =  sourcemodel.pos(sourcemodel.inside,:);

%% VS to ROIs
sourcemodel3 = sourcemodel2;
sourcemodel3.inside = sourcemodel.inside;

sourcedataproj = []; vs_roi = []; idx_pos = []; roi_coor = [];
for c = 1:length(vs.trial)
    a = vs.trial{c};
    tmp = zeros(size(sourcemodel3.pos,1),size(a,2));
    tmp(source_active.inside==1,:)= a;
    for ii = 1:length(atlas.tissuelabel)
        disp([num2str(ii),'/',num2str(length(atlas.tissuelabel))]);
        disp(atlas.tissuelabel(ii));
        l = Label== ii;
        idx_pos{ii} = sourcemodel3.pos(Label== ii,:);% xyz positions in mni coordinates
        if c==1
            roi_coor = [roi_coor;mean(idx_pos{ii},1)];
        end
        vs_source_act_sel = tmp(l,:);
        if ~isempty(idx_pos{ii}) && (vs_source_act_sel(1)>0)
            vs_source_act_sel(vs_source_act_sel(:,1)==0,:)=[];
            [u,s,v] = svd(vs_source_act_sel, 'econ');
            sourcedataproj(ii,:) = u(:,1)' * vs_source_act_sel;
        else
            sourcedataproj(ii,:) = zeros(1,size(vs.trial{1},2));
        end
        vs_roi{c} = sourcedataproj;
    end
end
vs_roi1 = [];
vs_roi1.trial = vs_roi;
vs_roi1.label = atlas.tissuelabel';
vs_roi1.time = vs.time;
vs_roi1.coor = roi_coor;

% figure; hold on;
% ft_plot_vol(vol, 'edgecolor', 'none', 'facealpha', 0.4);
% % ft_plot_mesh(sourcemodel2.pos(sourcemodel2.inside,:));
% ft_plot_mesh(vs_roi1.coor);
% ft_plot_mesh(source_active_temp.pos);

%% FFT comparision, VS vs. VS_ROI
% pause, close all,

% cfg = [];
% cfg.savefile = [];
% cfg.saveflag = 2;
% cfg.foilim = [2 50];
% cfg.plotflag  = 1;
% cfg.tapsmofrq       = 1;
% cfg.taper    = 'hanning';
% do_fft(cfg, vs); title('avg FFT of virtual sensors, voxels')
%
% cfg = [];
% cfg.savefile = [];
% cfg.saveflag = 2;
% cfg.foilim = [2 50];
% cfg.plotflag  = 1;
% cfg.tapsmofrq       = 1;
% cfg.taper    = 'hanning';
% do_fft(cfg, vs_roi1); title('avg FFT of virtual sensors, ROIs')