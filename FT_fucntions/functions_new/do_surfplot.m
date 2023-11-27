function do_surfplot(cfg)

% Author: vyoussofzadeh
% update: 04/25/22


% figure,
surf = cfg.surf;
d_in = cfg.d_in;

set(gcf, 'Position', cfg.position);
views = cfg.view;

for i=1:size(views,1)
    subplot(1,length(views),i)
    ft_plot_mesh(surf, 'vertexcolor', d_in);
%     camlight(80,-10);
%     camlight(-80,-10);
    colormap(cfg.color)
    colorbar off
    axis tight
    alpha(cfg.alpha)
    view(cfg.view(i,:));
    if ~isempty(cfg.coor)
        hold on
        plot3(cfg.coor(1),cfg.coor(2),cfg.coor(3),'r*')
    end
    hold off
end
title(cfg.title)

end