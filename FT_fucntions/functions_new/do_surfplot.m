function do_surfplot(cfg)

% Author: vyoussofzadeh
% update: 04/25/22

surf = cfg.surf;
d_in = cfg.d_in;

% figure('Position', cfg.position);  % Make sure a new figure is created
set(gcf, 'Position', cfg.position);
views = cfg.view;

% % %%
% % Number of colors in the colormap
% ncolors = 256;
% 
% % Split the colormap into two halves: one for negative curvature, one for positive
% half = ncolors / 2;
% 
% % Black to mid-gray (for negative curvature)
% neg_part = [ linspace(0, 0.5, half)' , ...  % R
%     linspace(0, 0.5, half)' , ...  % G
%     linspace(0, 0.5, half)' ];     % B
% 
% % Mid-gray to white (for positive curvature)
% pos_part = [ linspace(0.5, 1, half)' , ...  % R
%     linspace(0.5, 1, half)' , ...  % G
%     linspace(0.5, 1, half)' ];     % B
% 
% % Combine into a single colormap
% freesurferColormap = [neg_part; pos_part];
% 
% % Example usage in a plot:
% %  figure;
% %  trisurf(faces + 1, vertices(:,1), vertices(:,2), vertices(:,3), curvature_values, ...
% %      'EdgeColor', 'none', 'FaceAlpha', 1);
% %  colormap(freesurferColormap);
% %  caxis([-1 1]);         % Adjust depending on the range of your curvature data
%  axis equal; view(3); rotate3d on;


for i = 1:size(views,1)
    subplot(1,length(views),i)
    
    % Plot the mesh (FieldTrip version)
    ft_plot_mesh(surf, 'vertexcolor', d_in);
    
    % Turn off axis and set aspect ratio
    axis off
    axis tight
    
    % Set the colormap and colorbar options
    colormap(cfg.color)
    colorbar off
    
    % Adjust alpha (transparency)
    alpha(cfg.alpha)
    
    % Set the viewpoint
    view(views(i,:));
    
    % If coordinates are provided, plot them
    if ~isempty(cfg.coor)
        hold on
        plot3(cfg.coor(1), cfg.coor(2), cfg.coor(3), 'r*')
        hold off
    end
    
    
%     colormap(freesurferColormap);  % Use the custom color map
    
    axis off
    caxis([-1 1]);  % Adjust to fit the curvature range
    axis equal;
%     lighting phong;
%     camlight headlight;
    
end

% Add a title (it will appear over the last subplot)
sgtitle(cfg.title);

end


