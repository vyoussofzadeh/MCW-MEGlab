function do_surfplot(cfg,surf,source)

% Author: vyoussofzadeh
% update: 04/25/22

% figure,
set(gcf, 'Position', cfg.position);
views = cfg.view;

for i=1:length(views)
    subplot(1,length(views),i)
    handles = ft_plot_mesh(surf, 'vertexcolor', source);
    colormap(cfg.color)
    colorbar off
    axis tight
    view(cfg.view(i,:));
end
end