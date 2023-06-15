function do_conn_cluster(cfg_main, cln_data)

individual_grid = cfg_main.anat.individual_grid;
individual_headmodel = cfg_main.anat.individual_headmodel;
template_grid = cfg_main.anat.template_grid;
atlas = cfg_main.anat.atlas;
cov_matrix = cfg_main.cov_matrix;
seed = cfg_main.seed;
foi = cfg_main.foi;
pflag = cfg_main.pflag;

%% Virtual sensor
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
% template_grid.coordsys = 'mni';
% cfg = [];
% cfg.atlas      = atlas;
% cfg.roi        = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
% cfg.inputcoord = 'mni';
% mask           = ft_volumelookup(cfg, template_grid);

%%
% Label = zeros(length(individual_grid.pos),1);
% 
% for jj=1:length(atlas.tissuelabel)
%     
%     cfg = [];
%     cfg.atlas      = atlas;
%     cfg.roi        = atlas.tissuelabel{jj};  % here you can also specify a single label, i.e. single ROI
%     cfg.inputcoord = 'mni';
%     mask1           = ft_volumelookup(cfg, template_grid);
%     
%     individual_grid3 = individual_grid;
%     individual_grid3.inside = false(individual_grid3.dim);
%     individual_grid3.inside(mask1==1) = true;
%     
%     %                 hold on;
%     %                 ft_plot_mesh(individual_grid3.pos(individual_grid3.inside,:),'vertexcolor', [rand(1,1),rand(1,1),rand(1,1)]);
%     %             'edgecolor', [rand(1,1),rand(1,1),rand(1,1)]);
%     grid2 = individual_grid3;
%     grid2.inside = reshape(grid2.inside,length(grid2.pos),1);
%     
%     find(grid2.inside==1)
%     Label(grid2.inside==1,:)= jj;
%     
% end
% % figure, bar(Label(Label~=0))
% 
% %% Masking using aal atlas
% individual_grid2 = individual_grid;
% individual_grid2.inside = false(individual_grid2.dim);
% individual_grid2.inside(mask==1) = true;
% 
% if pflag.grid == 1
%     figure; hold on;
%     ft_plot_mesh(individual_headmodel.bnd, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.4;
%     ft_plot_mesh(individual_grid2.pos(individual_grid2.inside,:));
%     grid = individual_grid2;
%     grid.inside = reshape(grid.inside,length(grid.pos),1);
%     plot3(seed(1), seed(2), seed(3), 'm.','MarkerSize',80); % Sub2
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
% individual_grid3 = individual_grid2;
% individual_grid3.inside = template_grid.inside;
% 
% clc
% ft_progress('init', 'text',     'please wait ...');
% sourcedataproj = []; vs_roi = []; idx_pos = []; roi_coor = [];
% for c = 1:length(vs.trial)
%     ft_progress(c/length(vs.trial), 'extracting VS for trial %d from %d', c, length(vs.trial));
%     a = vs.trial{c};
%     tmp = zeros(size(individual_grid3.pos,1),size(a,2));
%     tmp(individual_grid3.inside==1,:)= a;
%     for ii = 1:length(atlas.tissuelabel)
% %         disp([num2str(ii),'/',num2str(length(atlas.tissuelabel))]);
% %         disp(atlas.tissuelabel(ii));
%         l = Label== ii;
%         idx_pos{ii} = individual_grid3.pos(Label== ii,:);% xyz positions in mni coordinates
%         if c==1
%             roi_coor = [roi_coor;mean(idx_pos{ii},1)];
%         end
%         vs_source_act_sel = tmp(l,:);
%         if ~isempty(idx_pos{ii}) && (vs_source_act_sel(1)>0)
%             vs_source_act_sel(vs_source_act_sel(:,1)==0,:)=[];
%             [u,~,~] = svd(vs_source_act_sel, 'econ');
%             sourcedataproj(ii,:) = u(:,1)' * vs_source_act_sel;
%         else
%             sourcedataproj(ii,:) = zeros(1,size(vs.trial{1},2));
%         end
%         vs_roi{c} = sourcedataproj;
%     end
%     pause(0.1);
% end
% ft_progress('close');
% 
% vs_roi1 = [];
% vs_roi1.trial = vs_roi;
% vs_roi1.label = atlas.tissuelabel';
% vs_roi1.time = vs.time;
% vs_roi1.coor = roi_coor;

%%
% seed_idx = knnsearch(vs_roi1.coor, seed);
% if pflag.grid_seed == 1    
%     figure; hold on;
%     ft_plot_mesh(individual_headmodel.bnd, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.4;
%     ft_plot_mesh(vs_roi1.coor);
%     hold on, plot3(seed(1), seed(2), seed(3), 'm.','MarkerSize',20); 
%     plot3(vs_roi1.coor(seed_idx,1), vs_roi1.coor(seed_idx,2), vs_roi1.coor(seed_idx,3), 'm.','MarkerSize',80, 'color', 'b'); % Sub2
%     
% end
% vs_roi1.label{seed_idx}

%%
grid_sel = individual_grid.pos(individual_grid.inside,:);
% seed_voxel_idx = knnsearch(grid_sel, seed);
[seed_idx, dist] = knnsearch(grid_sel, seed,'dist','cityblock','k',5);
if pflag.grid_seed == 1
    
    figure; hold on;
    ft_plot_mesh(individual_headmodel.bnd, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.4;
    ft_plot_mesh(individual_grid2.pos(individual_grid2.inside,:));
    grid = individual_grid2;
    grid.inside = reshape(grid.inside,length(grid.pos),1);
    for j=1:length(seed_idx)
        plot3(grid_sel(seed_idx(j),1), grid_sel(seed_idx(j),2), grid_sel(seed_idx(j),3), 'm.','MarkerSize',20, 'color', 'b'); % Sub2
    end
end

%%
for i=1:length(vs.trial)
   vs_all(i,:,:) = vs.trial{i};
   km(i,:) = kmeans(vs.trial{i}, 17);
end

% kk = kmeans(vs.trial{1}, 17);



vs_all1 = reshape(vs_all, [size(vs_all,2), size(vs_all,1)*size(vs_all,3)]);
km1 = kmeans(vs_all1, 17);

%%
conn_all = [];
for i=1:length(vs.trial)
   conn_all(i,:) = mean(do_conn(vs.trial{i}));
end

%%
reshapedData = reshape(vs_all, [], size(vs_all, 2));
standardizedData = zscore(reshapedData');
numComponents = 1; % Specify the desired number of principal components
[coefficients, scores, ~, ~, explained] = pca(standardizedData);
selectedComponents = coefficients(:, 1:numComponents);
reducedData = standardizedData * selectedComponents;


k = 5; % Number of clusters
[idx, centroids] = kmeans(reshapedData, k);

%%
net1 = [];
net1.pos = template_grid.pos;
net1.inside = template_grid.inside;
net1.clust = zeros(size(template_grid.pos,1),1);
net1.clust(net1.inside) = mean(km,1);
% net1.clust(net1.inside) = centroids(1,:);
% net.(netpar) = net.(netpar)./max(net.(netpar));
net1.clust(net1.inside) = mode(km);
% net1.clust(net1.inside) = reducedData;
% net1.clust(net1.inside) = mean(conn_all);

% net1.clust(net1.inside) = km1;

% mkm = mode(km);

close all
cfg = [];
cfg.subj = [num2str(foi), ' Hz'];
cfg.mask = 'clust';
cfg.thre = 0;
cfg.savepath = 'groupave';
cfg.colorbar = 2;
cfg.saveflag = 0;
cfg.colormap = colormap(flipud(pink));
cfg.colormap = colormap((jet));
% cfg.colormap = [];
cfg.surfinflated   = 'surface_inflated_both.mat';
% cfg.views = [90 0;0,90;-90 0;-180,-90];
cfg.views = [-90,0; 90,0];
% cfg.views = [-90,0;];
cfg.tit = '';
do_mapvis(cfg, net1);

% %% Conn analysis
% cfg            = [];
% cfg.output     = 'fourier';
% cfg.method     = 'mtmfft';
% cfg.foilim     = [1 30];
% cfg.tapsmofrq  = 2;
% cfg.keeptrials = 'yes';
% cfg.pad = 4;
% % freq    = ft_freqanalysis(cfg, vs_roi1);
% freq    = ft_freqanalysis(cfg, vs);
% 
% cfg         = [];
% cfg.method    = 'wpli_debiased';par = 'wpli_debiasedspctrm';
% source_conn = ft_connectivityanalysis(cfg, freq);
% 
% % cfg         = [];
% % cfg.method    = 'amplcorr'; par = 'amplcorrspctrm';
% % source_conn = ft_connectivityanalysis(cfg, freq);
% 
% % cfg         = [];
% % cfg.method  ='coh';
% % cfg.complex = 'absimag';
% % par = 'cohspctrm';
% % source_conn = ft_connectivityanalysis(cfg, freq);
% 
% %%
% source_conn1 = source_conn;
% % nonIdenticalIndices = 1:116 ~= seed_idx;
% % source_conn1.wpli_debiasedspctrm(:,nonIdenticalIndices,:) = 0;
% %
% % source_conn1.wpli_debiasedspctrm = source_conn1.wpli_debiasedspctrm(idx,:,:);
% for i=1:size(source_conn1.(par),3)
%     tmp = source_conn1.(par);
%     tmp(:,:,i) = tril(squeeze(tmp(:,:,i)));
% end
% source_conn1.(par) = tmp;
% 
% n = floor(size(source_conn1.(par),1)/15);
% % k=1;
% % y = round(linspace(1,size(source_conn1.(par),1),n));
% y =  seed_idx;
% 
% 
% if pflag.allconn == 1
%     for j=1:1%length(y)-1
%         figure
%         cfg           = [];
%         cfg.parameter = par;
%         cfg.zlim      = [0 1];
%         cfg.channel = y(j)-5:y(j)+5;%y(j):y(j+1);
%         ft_connectivityplot(cfg, source_conn1);
%         set(gcf, 'Position', [800   400   900   900]);
%         %             pause,
%     end
% end
% 
% %% network analysis
% netpar = 'eigenvector_cent'; net_thre = 0.5;
% % netpar = 'degrees'; net_thre = 0.5;
% 
% source_conn2  = source_conn1;
% % source_conn2.wpli_debiasedspctrm = source_conn2.wpli_debiasedspctrm(seed_idx,:,:);
% 
% cfg           = [];
% cfg.method    = netpar;
% cfg.parameter = par;
% cfg.threshold = net_thre;
% net = ft_networkanalysis(cfg,source_conn2);
% 
% f_range = round(net.freq);
% f_idx = 1:2:length(net.freq);
% % figure, imagesc(net.(netpar)),
% % set(gca,'Xtick', f_idx,'XtickLabel',f_range(f_idx));
% % set(gca,'Ytick', 1:size(net.(netpar),1),'YtickLabel',1:size(net.(netpar),1));
% % set(gca,'FontSize',9,'XTickLabelRotation',90);
% % set(gcf, 'Position', [800   400   500   1300]);
% 
% % for i=1:length(vs_roi1.label)
% %     roi_label{i} = [num2str(i), ':', vs_roi1.label{i}];
% % end
% % disp(roi_label');
% 
% % figure, plot(mean(net.(netpar))),
% figure, plot((net.(netpar)(seed_idx(1),:))),
% set(gca,'Xtick', f_idx,'XtickLabel',f_range(f_idx));
% set(gca,'Ytick', 1:size(net.(netpar),1),'YtickLabel',1:size(net.(netpar),1));
% set(gca,'FontSize',9,'XTickLabelRotation',90);
% % newString = strrep(vs_roi1.label{seed_idx}, '_', '-');
% % title(['seed', newString])
% set(gcf, 'Position', [1500   500   800  300]);
% 
% %%
% % net1 = net;
% % net1.degrees = net1.degrees(:,1);
% % net1.dimord = 'pos_pos';
% % net1.pos = vs_roi1.coor;
% % s_in = net1; mask = 'degrees'; thre = 0.6;
% % saveID = 'test';
% % Run_surfacemap_template
% 
% %% Plot AAL atlas
% % cfg = [];
% % cfg.funparameter = 'tissue';
% % cfg.method = 'surface';
% % cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% % cfg.projmethod     = 'nearest';
% % ft_sourceplot(cfg, atlas);
% % % view([0,90])
% % view([-90,0])
% 
% if pflag.aal ==1
%     cfg = [];
%     cfg.subj = ('Source map, mean');
%     cfg.mask = 'tissue';
%     cfg.thre = 0;
%     cfg.savepath = 'groupave';
%     cfg.colorbar = 2;
%     cfg.saveflag = 0;
%     cfg.colormap = [];
%     cfg.surfinflated   = 'surface_inflated_both.mat';
%     % cfg.views = [90 0;0,90;-90 0;-180,-90];
%     cfg.views = [-90,0; 90,0];
%     cfg.tit = '';
%     do_mapvis(cfg, atlas);
% end
% 
% %%
% disp(['from ', num2str(source_conn2.freq(1)), ' to ', num2str(source_conn2.freq(end)), 'Hz'])
% foi = input('frequncy range: ');
% 
% [~, idx_mn] = min(abs(foi(1) - source_conn2.freq));
% [~, idx_mx] = min(abs(foi(2) - source_conn2.freq));
% 
% % aedge =  squeeze(source_conn1.(par)(:,:,1));
% aedge =  mean(source_conn2.(par)(:,:,idx_mn:idx_mx),3);
% % aedge(nonIdenticalIndices,:) = 0;
% % size(aedge)
% aedge1 = zeros(size(aedge));
% aedge1(seed_idx,:) = aedge(seed_idx,:);
% aedge1(:,seed_idx) = aedge(:,seed_idx);
% 
% % ROI  = vs_roi1.coor;
% ROI = grid_sel;
% 
% aedge1(isnan(aedge1))=0; 
% tedge = (aedge1.* double(aedge1 > 0.1.*max(aedge1(:))));
% % figure, imagesc(tedge), colorbar
% 
% Verticies = individual_headmodel.bnd.pos;
% figure,
% % ft_plot_vol(individual_headmodel, 'facecolor', [0,0,0], 'edgecolor', 'none');
% % alpha 0.1;
% % view ([-196 56])
% plot3(Verticies(:,1),Verticies(:,2),Verticies(:,3),'color',[0.7,0.7,0.7]); box off, set(gca,'color','none'); axis off, axis image, rotate3d on
% % ft_plot_mesh(Verticies);
% hold on, h = plot3(ROI(:,1),ROI(:,2),ROI(:,3),'.g');
% set(h(1),'MarkerEdgeColor','g','MarkerFaceColor','g')
% hold on
% % k = 1;
% for i = 1:length(tedge)
%     for j = 1:length(tedge)
%         if tedge(i,j)> 0
%             p1 = [ROI(j,1),ROI(j,2),ROI(j,3)];
%             p2 = [ROI(i,1),ROI(i,2), ROI(i,3)];
%             pts = [p1; p2];
%             line(pts(:,1), pts(:,2), pts(:,3) ,'color','k');
%         end
%     end
% end
% view([0, 90])
% set(gcf, 'Position', [1500   500   300  300]);
% 
% %
% for i=1:length(atlas.tissuelabel)
%     roiid{i} = [num2str(i),': ',atlas.tissuelabel{i}];
% end
% disp(roiid')
% 
% %% color-coded network on aal atlas
% % close all
% % figure, imagesc(aedge);title(par); colorbar; colormap(flipud(pink));
% % set(gcf, 'Position', [1500   500   300  300]);
% 
% l = length(aedge);
% 
% % figure, bar(mean(aedge'),0.4); title('inflow-outflow');
% % set(gca,'Xtick', 1:l,'XtickLabel',1:l); box off, set(gca,'color','none');
% % set(gca,'FontSize',10,'XTickLabelRotation',90);
% % set(gcf, 'Position', [500   500   1500  300]);
% 
% % [a,b]  = sort(mean(aedge1),'descend');
% 
% % atlas1 = atlas;
% % for i = 1:l
% %     idx = atlas.tissue == b(i);
% %     atlas1.tissue(idx) = l-i;
% % end
% % unique(atlas.tissue);
% 
% source_conn3 = source_conn2;
% % source_conn3.(par) = mean(source_conn3.(par)(:,:,idx_mn:idx_mx),3);
% aedge1 = zeros(size(aedge));
% aedge1(seed_idx,:) = aedge(seed_idx,:);
% aedge1(:,seed_idx) = aedge(:,seed_idx);
% source_conn3.(par) = aedge1;
% 
% cfg           = [];
% cfg.method    = netpar;
% cfg.parameter = par;
% cfg.threshold = 0;
% net = ft_networkanalysis(cfg,source_conn3);
% net.pos     = template_grid.pos;
% net.dim     = template_grid.dim;
% net.inside  = template_grid.inside;
% tmp = net.(netpar);
% % tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:)));
% net.(netpar) = zeros(size(net.pos,1),1);
% net.(netpar)(net.inside) = tmp';
% 
% close all
% cfg = [];
% cfg.subj = [num2str(foi), ' Hz'];
% cfg.mask = netpar;
% cfg.thre = 0;
% cfg.savepath = 'groupave';
% cfg.colorbar = 2;
% cfg.saveflag = 0;
% % cfg.colormap = colormap(flipud(pink));
% cfg.colormap = [];
% cfg.surfinflated   = 'surface_inflated_both.mat';
% % cfg.views = [90 0;0,90;-90 0;-180,-90];
% cfg.views = [-90,0; 90,0];
% % cfg.views = [-90,0;];
% cfg.tit = '';
% do_mapvis(cfg, net);
% % colorbar
% 
% % roiid(b(1:10))'
% 
% %% color-coded seed on aal atlas
% seed_map = atlas;
% seed_map.tissue = zeros(size(atlas.tissue));
% for i = 1:l
%     idx = atlas.tissue == seed_idx;
%     seed_map.tissue(idx) = 1;
% end
% 
% cfg = [];
% cfg.subj = ('Source map, mean');
% cfg.mask = 'tissue';
% cfg.thre = 0;
% cfg.savepath = 'groupave';
% cfg.colorbar = 2;
% cfg.saveflag = 0;
% cfg.colormap = [];
% cfg.surfinflated   = 'surface_inflated_both.mat';
% % cfg.views = [90 0;0,90;-90 0;-180,-90];
% cfg.views = [-90,0; 90,0];
% cfg.tit = '';
% do_mapvis(cfg, seed_map);
% 
% %%
% out = [];
% out.network = net;
% out.foi = foi;
% out.source_active = source_active;
% out.seed_idx = seed_idx;
% out.seed_roi = vs_roi1.label{seed_idx};
% out.network_atlas = atlas1;
% out.source_conn = source_conn2;