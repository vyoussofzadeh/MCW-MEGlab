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
        print(gcf,fullfile(cfg_main.savepath, [cfg_main.subj,'_',num2str(j)]),'-dpng');
        %         saveas(gcf,[cfg_main.subj,'_',num2str(j)]); % matlab fig
        %         end
    end
    %     pause(1)
    rotate3d on
    set(gcf, 'Position', [1500   500   300  300]);
end

if cfg_main.saveflag ==1
    %% Combine images
    a = [];
    cd(cfg_main.savepath)
    for j=1:size(cfg_main.views,1)
        b=imread(fullfile(cfg_main.savepath, [cfg_main.subj,'_',num2str(j),'.png']));
        if isa(b,'uint8'), b=double(b)/255; end
        if max(b(:))>1, b=double(b)/double(max(b(:))); end
        a{j}=double(b);
    end
    
    ncut = 16;
    domosaiccrop = 1;
    if domosaiccrop
        cropt_idx={};
        for n=1:numel(a)
            cropt=any(any(diff(a{n},1,2),2),3);
            cropt_idx{n,1}=max(1,sum(~cumsum(cropt))-ncut):size(a{n},1)-max(0,sum(~cumsum(flipud(cropt)))-ncut);
            cropt=any(any(diff(a{n},1,1),1),3);
            cropt_idx{n,2}=max(1,sum(~cumsum(cropt))-ncut):size(a{n},2)-max(0,sum(~cumsum(flipud(cropt)))-ncut);
        end
    end
    
    if domosaiccrop
        cropt_idx1234=cropt_idx{1,1}; for n=2:numel(a), cropt_idx1234=union(cropt_idx1234,cropt_idx{n,1}); end
        ta=[]; for n=1:numel(a), ta=cat(2,ta,a{n}(cropt_idx1234,cropt_idx{n,2},:)); end; a=ta;
        %cropt_idx1234=union(union(union(cropt_idx{1,1},cropt_idx{2,1}),cropt_idx{3,1}),cropt_idx{4,1});
        %a=[a{1}(cropt_idx1234,cropt_idx{1,2},:),a{2}(cropt_idx1234,cropt_idx{2,2},:),a{3}(cropt_idx1234,cropt_idx{3,2},:),a{4}(cropt_idx1234,cropt_idx{4,2},:)];
    else
        a=cat(2,a{:});
        %a=[a{1},a{2},a{3},a{4}];
    end
    
    imwrite(a,fullfile(cfg_main.savepath, [cfg_main.subj,'.png']));
    
    %%
    
    for i=1:j
        delete(fullfile(cfg_main.savepath, [cfg_main.subj,'_',num2str(i),'.png']))
    end
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

% figure,
% subplot(1,2,1)
% copyobj(allchild(get(1,'CurrentAxes')),gca); axis off, axis square
% view(cfg_main.views);
% camlight;
% material dull;
% colorbar off
% axis tight;
% subplot(1,2,2)
% copyobj(allchild(get(1,'CurrentAxes')),gca); axis off, axis square
% view([90,0]);
% % set(gcf, 'Position', [1500   500   300  300]);
% camlight;
% material dull;
% colorbar off
% % set(gcf,'Color','k')
% axis tight;
% set(gcf, 'Position', [1500   500   700  200]);
% tightfig
%
% ax = subplot(1,2,2);
% pos = get(ax, 'position');
% pos(1) = 0.1;  % set the x-position to 0.1
% set(ax, 'position', pos);

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
