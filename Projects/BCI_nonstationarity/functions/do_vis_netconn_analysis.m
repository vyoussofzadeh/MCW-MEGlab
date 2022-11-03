function do_vis_netconn_analysis(mcfg, source_conn1)

disp(['from ', num2str(source_conn1.freq(1)), ' to ', num2str(source_conn1.freq(end)), 'Hz'])
foi = input('frequncy range e.g, [5,25]: ');

[mn, idx_mn] = min(abs(foi(1) - source_conn1.freq));
[mx, idx_mx] = min(abs(foi(2) - source_conn1.freq));

% aedge =  squeeze(source_conn1.(par)(:,:,1));
aedge =  mean(source_conn1.(mcfg.par)(:,:,idx_mn:idx_mx),3);
% size(aedge)

ROI  = mcfg.coor;

aedge(isnan(aedge)) = 0; tedge = (aedge.* double(aedge > 0.8.*max(aedge(:))));
% figure, imagesc(tedge), colorbar

Verticies = mcfg.individual_headmodel.bnd.pos;
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
set(gcf, 'Position', [1000   500   300  300]);

if mcfg.flag.network_exportfig == 1
    print(['Conn_', num2str(foi(1)),'-', num2str(foi(2)), ' Hz', ],'-dpng');
end

%%
atlas = mcfg.atlas;
for i=1:length(atlas.tissuelabel)
    roiid{i} = [num2str(i),': ',atlas.tissuelabel{i}];
end
disp(roiid')

%%
l = length(aedge);

% figure, 
% bar(mean(aedge),0.4); title('mean conn values');
% gg = gcf;
% set(gca,'Xtick', 1:l,'XtickLabel',1:l); box off, set(gca,'color','none');
% set(gca,'FontSize',10,'XTickLabelRotation',90);
% set(gg, 'Position', [500   500   1500  300]);

[a,b]  = sort(mean(aedge),'descend');
% b(b < 0.7.*max(b)) = 0;

atlas1 = atlas;
for i = 1:15%l
    idx = find(atlas.tissue == b(i));
    atlas1.tissue(idx) = l-i;
end
unique(atlas.tissue);

cfg = [];
cfg.subj = [num2str(foi(1)),'-', num2str(foi(2)), ' Hz', ];
cfg.mask = 'tissue';
cfg.thre = 0;
cfg.savepath = 'groupave';
cfg.colorbar = 2;
cfg.saveflag = 2;
% cfg.colormap = colormap(flipud(pink));
% cfg.colormap = [];
cfg.colormap = flipud(brewermap(64,'RdBu'));
cfg.surfinflated   = 'surface_inflated_both.mat';
% cfg.views = [90 0;0,90;-90 0;-180,-90];
cfg.views = [-90,0; 90,0];
do_mapvis(cfg, atlas1);
% colorbar

roiid(b(1:10))'