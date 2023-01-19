%% Virtual sensor
% pause,
% close all,

cfg = [];
cfg.method = 'lcmv';
cfg.grid = individual_grid;
cfg.headmodel = individual_headmodel;
cfg.lcmv.lambda = '10%';
cfg.lcmv.keepfilter = 'yes';
source_whole = ft_sourceanalysis(cfg, cov_matrix);

cfg =[];
cfg.method = 'lcmv';
cfg.grid = individual_grid;
cfg.grid.filter = source_whole.avg.filter;
cfg.headmodel = individual_headmodel;
cfg.rawtrial = 'yes';
cfg.keeptrials = 'yes';
cfg.lcmv.projectmom = 'yes';
cfg.lcmv.fixedori = 'yes';
source_active = ft_sourceanalysis(cfg, cln_data);

%%
source_whole_temp = source_whole;
source_whole_temp.pos     = template_grid.pos;
source_whole_temp.dim     = template_grid.dim;
source_whole_temp.inside  = template_grid.inside;

source_active_temp = source_active;
source_active_temp.pos     = template_grid.pos;
source_active_temp.dim     = template_grid.dim;
source_active_temp.inside  = template_grid.inside;

%%
template_grid.coordsys = 'mni';
cfg = [];
cfg.atlas      = atlas;
cfg.roi        = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
cfg.inputcoord = 'mni';
mask           = ft_volumelookup(cfg, template_grid);

%%
% pause,
% close all,
%         figure;
%         ft_plot_vol(individual_headmodel, 'edgecolor', 'none', 'facealpha', 0.4);
Label = zeros(length(individual_grid.pos),1);

for jj=1:length(atlas.tissuelabel)
    
    cfg = [];
    cfg.atlas      = atlas;
    cfg.roi        = atlas.tissuelabel{jj};  % here you can also specify a single label, i.e. single ROI
    cfg.inputcoord = 'mni';
    mask1           = ft_volumelookup(cfg, template_grid);
    
    
    
    individual_grid3 = individual_grid;
    individual_grid3.inside = false(individual_grid3.dim);
    individual_grid3.inside(mask1==1) = true;
    
    %                 hold on;
    %                 ft_plot_mesh(individual_grid3.pos(individual_grid3.inside,:),'vertexcolor', [rand(1,1),rand(1,1),rand(1,1)]);
    %             'edgecolor', [rand(1,1),rand(1,1),rand(1,1)]);
    grid2 = individual_grid3;
    grid2.inside = reshape(grid2.inside,length(grid2.pos),1);
    
    find(grid2.inside==1)
    Label(grid2.inside==1,:)= jj;
    
end
figure, bar(Label(Label~=0))

%% Masking using aal atlas
% pause, close all,
individual_grid2 = individual_grid;
individual_grid2.inside = false(individual_grid2.dim);
individual_grid2.inside(mask==1) = true;

figure; hold on;
ft_plot_vol(individual_headmodel, 'edgecolor', 'none', 'facealpha', 0.4);
ft_plot_mesh(individual_grid2.pos(individual_grid2.inside,:));
grid = individual_grid2;
grid.inside = reshape(grid.inside,length(grid.pos),1);

%%
%         cfg = [];
%         cfg.method    = 'nearest';
%         cfg.parameter = 'tissue';
%         atlas2 = ft_sourceinterpolate(cfg,atlas,individual_grid);
%
%         %%
%         cfg = [];
%         cfg.method       = 'eig';
%         cfg.parcellation = 'tissue';
%         ft_sourceparcellate(cfg,source_active,atlas2)

%% VS extraction
% % pause, close all,
% 
% active_nodes = find(individual_grid.inside==1);
% % create the trial field from source anal output
% vs = [];
% node = [];
% for i = 1:length(source_active.trial) % i = number trials
%     for x = 1:length(active_nodes) % x = number nodes
%         node(x,:) = source_active.trial(i).mom{active_nodes(x)};
%     end
%     vs.trial{i} = node;
% end
% 
% %
% % set up labels
% label_temp = num2str(active_nodes);
% label_temp = cellstr(label_temp);
% vs.label = label_temp;
% 
% % set up time component
% for i = 1:length(source_active.trial)
%     vs.time{1, i} = source_active.time;
% end

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

%
% set up labels
label_temp = num2str(active_nodes);
label_temp = cellstr(label_temp);
vs.label = label_temp;

% set up time component
for i = 1:length(source_active_temp.trial)
    vs.time{1, i} = source_active_temp.time;
end

%% VS to ROIs
% sourcedataproj = []; vs_roi = []; idx_pos = []; roi_coor = [];
% for c = 1:length(vs.trial)
%     a = vs.trial{c};
%     tmp = zeros(size(template_grid.pos,1),size(a,2));
%     tmp(template_grid.inside==1,:)= a;
%     for ii = 1:length(atlas.tissuelabel)
%         disp([num2str(ii),'/',num2str(length(atlas.tissuelabel))]);
%         disp(atlas.tissuelabel(ii));
%         l = find (Label== ii);
%         idx_pos{ii} = template_grid.pos(Label== ii,:);% xyz positions in mni coordinates
%         if c==1
%             roi_coor = [roi_coor;mean(idx_pos{ii},1)];
%         end
%         vs_source_act_sel = tmp(l,:);
%         if ~isempty(idx_pos{ii}) && (vs_source_act_sel(1)>0)
%             vs_source_act_sel(vs_source_act_sel(:,1)==0,:)=[];
%             [u,s,v] = svd(vs_source_act_sel, 'econ');
%             sourcedataproj(ii,:) = u(:,1)' * vs_source_act_sel;
%         else
%             sourcedataproj(ii,:) = zeros(1,size(vs.trial{1},2));
%         end
%         vs_roi{c} = sourcedataproj;
%     end
% end

individual_grid3 = individual_grid2;
individual_grid3.inside = template_grid.inside;

sourcedataproj = []; vs_roi = []; idx_pos = []; roi_coor = [];
for c = 1:length(vs.trial)
    a = vs.trial{c};
    tmp = zeros(size(individual_grid3.pos,1),size(a,2));
    tmp(individual_grid3.inside==1,:)= a;
    for ii = 1:length(atlas.tissuelabel)
        disp([num2str(ii),'/',num2str(length(atlas.tissuelabel))]);
        disp(atlas.tissuelabel(ii));
        l = Label== ii;
        idx_pos{ii} = individual_grid3.pos(Label== ii,:);% xyz positions in mni coordinates
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


figure; hold on;
ft_plot_vol(individual_headmodel, 'edgecolor', 'none', 'facealpha', 0.4);
% ft_plot_mesh(individual_grid2.pos(individual_grid2.inside,:));
ft_plot_mesh(vs_roi1.coor);
% ft_plot_mesh(source_active_temp.pos);

%% FFT comparision, VS vs. VS_ROI
% pause, close all,

cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 50];
cfg.plotflag  = 1;
cfg.tapsmofrq       = 1;
cfg.taper    = 'hanning';
vy_fft(cfg, vs);

cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 50];
cfg.plotflag  = 1;
cfg.tapsmofrq       = 1;
cfg.taper    = 'hanning';
vy_fft(cfg, vs_roi1);

%% Conn analysis
% pause, close all,
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [1 30];
cfg.tapsmofrq  = 2;
cfg.keeptrials = 'yes';
cfg.pad = 4;
freq    = ft_freqanalysis(cfg, vs_roi1);

cfg         = [];
cfg.method    = 'wpli_debiased';
source_conn = ft_connectivityanalysis(cfg, freq);
par = 'wpli_debiasedspctrm';

% cfg         = [];
% cfg.method    = 'amplcorr';
% source_conn = ft_connectivityanalysis(cfg, freq);
% par = 'amplcorrspctrm';

%%
source_conn1 = source_conn;
for i=1:size(source_conn1.(par),3)
    tmp = source_conn1.(par);
    tmp(:,:,i) = tril(squeeze(tmp(:,:,i)));
end
source_conn1.(par) = tmp;

n = floor(size(source_conn1.(par),1)/15); 
% k=1;
y = round(linspace(1,size(source_conn1.(par),1),n));

for j=1:1%length(y)-1
    figure
    cfg           = [];
    cfg.parameter = par;
    cfg.zlim      = [0 1];
    cfg.channel = y(j):y(j+1);
    ft_connectivityplot(cfg, source_conn1);
    set(gcf, 'Position', [800   400   900   900]);
    %             pause,
end

%%
% mfreq = mean(freq.fourierspctrm,2)
% 
% grangercfg = [];
% grangercfg.method  = 'granger';
% grangercfg.granger.conditional = 'no';
% grangercfg.granger.sfmethod = 'bivariate';
% 
% gdata = [];
% gdata.g1_bivar_reg      = ft_connectivityanalysis(grangercfg, freq);
% 
% %%
% conn = gdata.g1_bivar_reg;
% mask = 'grangerspctrm';
% 
% % test = permute(conn.(mask),[116,116,113]);
% 
% % [r c d] = size(conn.(mask));       % Originally a 340x740x144 array
% % img = reshape(data,r*c,d);  % This gives me a 251600x144 array
% test = conn.(mask);
% Y = permute(reshape(test, 115, 116, 113), [1, 2, 3]);
% 
% % [r,c] = size(conn.(mask));
% % nlay  = 3;
% % out   = reshape(conn.(mask),[116,116,113]);
% 
% source_conn1 = Y;
% for i=1:size(source_conn1.(mask),3)
%     source_conn1.(mask)(:,:,i) = tril(squeeze(source_conn1.(mask)(:,:,i)));
% end
% 
% n = floor(size(source_conn1.(mask),1)/15); k=1;
% y = round(linspace(1,size(source_conn1.(mask),1),n));
% 
% for j=1:1%length(y)-1
%     figure
%     cfg           = [];
%     cfg.parameter = mask;
%     cfg.zlim      = [0 1];
%     cfg.channel = y(j):y(j+1);
%     ft_connectivityplot(cfg, source_conn1);
%     set(gcf, 'Position', [800   400   900   900]);
%     %             pause,
% end

%% network analysis
% pause, close all,
cfg           = [];
cfg.method    = 'degrees';
cfg.parameter = par;
cfg.threshold = .5;
net = ft_networkanalysis(cfg,source_conn1);

f_range = round(net.freq);
f_idx = 1:2:length(net.freq);
figure, imagesc(net.degrees),
set(gca,'Xtick', f_idx,'XtickLabel',f_range(f_idx));
set(gca,'Ytick', 1:size(net.degrees,1),'YtickLabel',1:size(net.degrees,1));
set(gca,'FontSize',9,'XTickLabelRotation',90);
set(gcf, 'Position', [800   400   500   1300]);

for i=1:length(vs_roi1.label)
    roi_label{i} = [num2str(i), ':', vs_roi1.label{i}];
end
disp(roi_label');

figure, plot(mean(net.degrees)),
set(gca,'Xtick', f_idx,'XtickLabel',f_range(f_idx));
set(gca,'Ytick', 1:size(net.degrees,1),'YtickLabel',1:size(net.degrees,1));
set(gca,'FontSize',9,'XTickLabelRotation',90);

%%
% net1 = net;
% net1.degrees = net1.degrees(:,1);
% net1.dimord = 'pos_pos';
% net1.pos = vs_roi1.coor;
% s_in = net1; mask = 'degrees'; thre = 0.6;
% saveID = 'test';
% Run_surfacemap_template

%% Plot AAL atlas
% cfg = [];
% cfg.funparameter = 'tissue';
% cfg.method = 'surface';
% cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% cfg.projmethod     = 'nearest'; 
% ft_sourceplot(cfg, atlas);
% % view([0,90])
% view([-90,0])

% cfg = [];
% cfg.subj = ('Source map, mean');
% cfg.mask = 'tissue';
% cfg.thre = 0;
% cfg.savepath = 'groupave';
% cfg.colorbar = 2;
% cfg.saveflag = 0;
% cfg.surfinflated   = 'surface_inflated_both.mat';
% % cfg.views = [90 0;0,90;-90 0;-180,-90];
% cfg.views = [-90,0];
% vy_mapvisualisation2(cfg, atlas);

%%
disp(['from ', num2str(source_conn1.freq(1)), ' to ', num2str(source_conn1.freq(end)), 'Hz'])
foi = input('frequncy range: ');

[mn, idx_mn] = min(abs(foi(1) - source_conn1.freq));
[mx, idx_mx] = min(abs(foi(2) - source_conn1.freq));

% aedge =  squeeze(source_conn1.(par)(:,:,1));
aedge =  mean(source_conn1.(par)(:,:,idx_mn:idx_mx),3);
% size(aedge)

ROI  = vs_roi1.coor;

aedge(isnan(aedge))=0; tedge = (aedge.* double(aedge > 0.5.*max(aedge(:))));
% figure, imagesc(tedge), colorbar

Verticies = individual_headmodel.bnd.pos;
figure,
% ft_plot_vol(individual_headmodel, 'facecolor', [0,0,0], 'edgecolor', 'none');
% alpha 0.1;
% view ([-196 56])
plot3(Verticies(:,1),Verticies(:,2),Verticies(:,3),'color',[0.7,0.7,0.7]); box off, set(gca,'color','none'); axis off, axis image, rotate3d on
hold on, h = plot3(ROI(:,1),ROI(:,2),ROI(:,3),'.g');
set(h(1),'MarkerEdgeColor','g','MarkerFaceColor','g')
hold on
k = 1;
for i = 1:length(tedge)
    for j = 1:length(tedge)
        if tedge(i,j)> 0
            p1 = [ROI(j,1),ROI(j,2),ROI(j,3)];
            p2 = [ROI(i,1),ROI(i,2), ROI(i,3)];
            pts = [p1; p2];
            line(pts(:,1), pts(:,2), pts(:,3) ,'color','k');
        end
    end
end
view([90 0])

%%
for i=1:length(atlas.tissuelabel)
    roiid{i} = [num2str(i),': ',atlas.tissuelabel{i}];
end
disp(roiid')

%%
% close all
% figure, imagesc(aedge);title(par); colorbar; colormap(flipud(pink));
% set(gcf, 'Position', [1500   500   300  300]);

l = length(aedge);

figure, bar(mean(aedge),0.4); title('inflow-outflow');
set(gca,'Xtick', 1:l,'XtickLabel',1:l); box off, set(gca,'color','none');
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gcf, 'Position', [500   500   1500  300]);

[a,b]  = sort(mean(aedge),'descend');

atlas1 = atlas;
for i = 1:l
    idx = find(atlas.tissue == b(i));
%     if i==1
        atlas1.tissue(idx) = l-i;
%     else
%         atlas1.tissue(idx) = 1;
%     end
end
unique(atlas.tissue);

cfg = [];
cfg.subj = [num2str(foi), ' Hz'];
cfg.mask = 'tissue';
cfg.thre = 0;
cfg.savepath = 'groupave';
cfg.colorbar = 2;
cfg.saveflag = 0;
% cfg.colormap = colormap(flipud(pink)); 
cfg.colormap = []; 
cfg.surfinflated   = 'surface_inflated_both.mat';
% cfg.views = [90 0;0,90;-90 0;-180,-90];
cfg.views = [-90,0; 90,0];
vy_mapvisualisation2(cfg, atlas1);
% colorbar

roiid(b(1:10))'

%%  Freq analysis and MVAR modelling (parametric)
% cfg         = [];
% cfg.order   = 20;
% % cfg.order   = aic;
% % cfg.toolbox = 'bsmart';
% Dmdata       = ft_mvaranalysis(cfg, vs_roi1);
% 
% % ft_hastoolbox('bsmart',1)
% 
% cfg        = [];
% cfg.method = 'mvar';
% Dmfreq      = ft_freqanalysis(cfg, Dmdata);
% 
% %%
% method = 'pdc'; par = 'pdcspctrm';
% cfg           = [];
% cfg.method    = method;
% par_conn           = ft_connectivityanalysis(cfg, Dmfreq);


%%
% vol1 = ft_convert_units(individual_headmodel, 'cm');
% Verticies = individual_headmodel.bnd.pos;
% figure,
% plot3(Verticies(:,1),Verticies(:,2),Verticies(:,3),'color',[0.7,0.7,0.7]);
% hold on
% h = plot3(ROI(:,1),ROI(:,2),ROI(:,3),'.g');
% % set(h(1),'MarkerEdgeColor',[1 0.48 0.30],'MarkerFaceColor','g')
% set(h(1),'MarkerEdgeColor','g','MarkerFaceColor','g')
% % set(h(2),'MarkerEdgeColor','none','MarkerFaceColor','g')
% box off
% set(gca,'color','none');
% axis off
% axis image
% rotate3d on
% hold on
% view([-90,90])
% % view ([180 90])

%% FOI_ connectvitiy (incomplete)
%         cfg = [];
%         cfg.channel = [78:79];
%         source_conn1 = ft_selectdata(cfg, source_conn1)
%
%         close all
%         cfg           = [];
%         cfg.method    = 'degrees';
%         cfg.parameter = '(par)';
%         cfg.threshold = .5;
%         cfg.channels = [78:79];
%         net1 = ft_networkanalysis(cfg,source_conn1);
%
%
%         source_conn2 = source_conn1;
%         source_conn2.(par) =
%
%         close all
%         cfg           = [];
%         cfg.method    = 'degrees';
%         cfg.parameter = '(par)';
%         cfg.threshold = .5;
%         net = ft_networkanalysis(cfg,source_conn1);
%
%         %
%         figure, plot(mean(net.degrees(78:79,:))),
%         set(gca,'Xtick', f_idx,'XtickLabel',f_range(f_idx));
%         set(gca,'Ytick', 1:size(net.degrees,1),'YtickLabel',1:size(net.degrees,1));
%         set(gca,'FontSize',9,'XTickLabelRotation',90);


%% Hypothesis testing
%         net.pos     = template_grid.pos;
%         net.dim     = template_grid.dim;
%         net.inside  = template_grid.inside;
%
%
%         cfg = [];
%         cfg.interpmethod = 'nearest';
%         cfg.parameter    = 'tissue';
%         parcel_atlas = ft_sourceinterpolate(cfg, atlas, net);
%
%         x = find(ismember(atlas.tissuelabel,'Heschl_L'));
%         indxHGL = find(parcel_atlas.tissue == x);
%
%         x=find(ismember(atlas.tissuelabel,'Heschl_R'));
%         indxHGR = find(parcel_atlas.tissue==x);
%
%         x=find(ismember(atlas.tissuelabel,'Cingulum_Mid_L'));
%         indxCML = find(parcel_atlas.tissue==x);




%         cfg = [];
%         cfg.channel = [13,15];% {'Frontal_Inf_Tri_L'}
%         cfg.avgoverchan = 'yes';
%         vs_1 = ft_selectdata(cfg,vs_roi1);
%
%         cfg = [];
%         cfg.channel = [10,12];% {'Frontal_Inf_Tri_L'}
%         cfg.avgoverchan = 'yes';
%         vs_2 = ft_selectdata(cfg,vs_roi1);
%
%         virtsensparcel=ft_appenddata([],vs_1,vs_2);

%%
%         cfg            = [];
%         cfg.output     = 'fourier';
%         cfg.method     = 'mtmfft';
%         cfg.foilim     = [10 30];
%         cfg.tapsmofrq  = 2;
%         cfg.keeptrials = 'yes';
%         cfg.pad = 4;
%         freq    = ft_freqanalysis(cfg, vs_roi1);
%