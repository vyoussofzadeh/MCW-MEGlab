function do_mapvis(cfg_main, input)

cfg                = [];
cfg.method         = 'ortho';
cfg.funparameter   = cfg_main.mask;
cfg.funcolorlim    = 'maxabs';
cfg.opacitymap     = 'rampup';
cfg.crosshair      = 'no';
cfg.camlight       = 'no';
cfg.funcolormap    =  cfg_main.colormap;
cfg.projthresh     = cfg_main.thre;

cfg.method = 'surface';
cfg.surfinflated   = cfg_main.surfinflated;
for j=1:size(cfg_main.views,1)
    
    ft_sourceplot(cfg, input);
    view(cfg_main.views(j,:));
    camlight;
    material dull;
    colorbar off
    set(gcf,'name',cfg_main.tit,'numbertitle','off')
    if cfg_main.saveflag ==1
        %         pause(1)
        set(gcf,'Color','k')
        set(gcf, 'Position', [1500   500   300  300]);
        %         if ~exist([cfg_main.subj,'_',num2str(j),'.png'],'file')
        print(gcf,[cfg_main.subj,'_',num2str(j)],'-dpng');
%         saveas(gcf,[cfg_main.subj,'_',num2str(j)]); % matlab fig
        %         end
    end
    %     pause(1)
    rotate3d on
    set(gcf, 'Position', [1500   500   300  300]);
end

%%
% Combine the two figures into a single figure with 1 row and 2 columns
% figure(j+1);
% subplot(1,2,1);
% copyobj(allchild(get(1,'CurrentAxes')),gca); axis off, axis square
% title('Combined Plot 1');
% subplot(1,2,2);
% copyobj(allchild(get(1,'CurrentAxes')),gca); axis off, axis square
% title('Combined Plot 2');
% 
% 
% %%
% figure(j+1);
% 
% for j=1:size(cfg_main.views,1)
%     
%     subplot(1,size(cfg_main.views,1),j);
%     copyobj(allchild(get(1,'CurrentAxes')),gca);
%     
%     axis tight;
%     axis square;
%     axis off,
%     
%     %     axis square
%     %     ft_sourceplot(cfg, input);
%     view(cfg_main.views(j,:));
%     %     camlight;
%     %     material dull;
%     %     colorbar off
%     %     pause(1)
%     %     ff = gcf;
%     %     set(gcf,'name',cfg_main.tit,'numbertitle','off')
%     %     if cfg_main.saveflag ==1
%     % %         pause(1)
%     %         set(gcf,'Color','k')
%     %         set(gcf, 'Position', [1500   500   300  300]);
%     % %         if ~exist([cfg_main.subj,'_',num2str(j),'.png'],'file')
%     %             print(gcf,[cfg_main.subj,'_',num2str(j)],'-dpng');
%     %             saveas(gcf,[cfg_main.subj,'_',num2str(j)]);
%     % %         end
%     %     end
%     % %     pause(1)
%     %     rotate3d on
%     %     set(gcf, 'Position', [1500   500   300  300]);
% end
% set(gcf, 'Position', [1500   500   1000  200]);
