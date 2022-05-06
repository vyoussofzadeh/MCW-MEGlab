function vy_mapvisualisation2(cfg_main, input)

% ,funparameter,thre,savepath

cfg                = [];
cfg.method         = 'ortho';
cfg.funparameter   = cfg_main.mask;
cfg.funcolorlim    = 'maxabs';
cfg.opacitymap     = 'rampup';
cfg.crosshair      = 'no';
cfg.camlight       = 'no';
% cfg.funcolormap    =  brewermap(256, '*RdYlBu');
cfg.funcolormap    =  cfg_main.colormap;
cfg.projthresh     = cfg_main.thre;

cfg.method = 'surface';
% cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% cfg.surfinflated   = 'surface_pial_both.mat';
cfg.surfinflated   = cfg_main.surfinflated;
% 'surface_inflated_both.mat';
% cfg.surfinflated   = 'surface_inflated_right.mat';
for j=1:size(cfg_main.views,1)
    
    ft_sourceplot(cfg, input);
    % view ([-70 20 50])
    % light ('Position',[-70 20 50]);
    view(cfg_main.views(j,:));
    camlight;
    material dull;
    colorbar off
    set(gcf,'name',cfg_main.subj,'numbertitle','off')
    if cfg_main.saveflag ==1
        %     hcp_write_figure([savepath{1},'.png'], gcf, 'resolution', 300);
        pause(0.5)
        %     print(cfg_main.savepath{1},'-depsc')
        %     print(cfg_main.savepath{1},'-dpng')
        %     print([cfg_main.subj,'_L'],'-dpdf');
        %     print([cfg_main.subj,'_L'],'-dpng');
%         set(gcf,'Color','w')
        set(gcf,'Color','k')
        set(gcf, 'Position', [1500   500   300  300]);
        if ~exist([cfg_main.subj,'_',num2str(j),'.png'],'file')
            print([cfg_main.subj,'_',num2str(j)],'-dpng');
            %     export_fig([cfg_main.subj,'.pdf'],'-append')
        end
    end
    rotate3d on
    set(gcf, 'Position', [1500   500   300  300]);
end

% if cfg_main.colorbar ==2
%     colorbar off
% end
% set(gcf,'name',cfg_main.subj,'numbertitle','off')
% title('test');



% % figure,
% % subplot 211
% ft_sourceplot(cfg, input);
% % view ([70 20 50])
% % light ('Position',[70 20 50])
% view([-90 0]); camlight; material dull;
% % if cfg_main.colorbar ==2
% %     colorbar off
% % end
% set(gcf,'name',cfg_main.subj,'numbertitle','off')
% colorbar off
% if cfg_main.saveflag== 1
%     %     hcp_write_figure([savepath{2},'.png'], gcf, 'resolution', 300);
%     pause(0.5)
%     %     print(cfg_main.savepath{2},'-depsc')
%     %     print(cfg_main.savepath{2},'-dpng')
%     %     print([cfg_main.subj,'_R'],'-dpdf');
%     set(gcf,'Color','w')
%     print([cfg_main.subj,'_L'],'-dpng');
%     %     export_fig([cfg_main.subj,'.pdf'],'-append')
% end



