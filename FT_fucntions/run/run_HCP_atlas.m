% load in LH
[~, lhlab,lhctable]=read_annotation('lh.HCP-MMP1.annot');
[lhvtx,lhfaces]=read_surf('lh.pial');
[~,lhix]=ismember(lhlab,lhctable.table(:,5));

% load in RH
[~, rhlab,rhctable]=read_annotation('rh.HCP-MMP1.annot');
[rhvtx,rhfaces]=read_surf('rh.pial');
[~,rhix]=ismember(rhlab,rhctable.table(:,5));

if flag.plotatlas ==1
    close all
    % visualize Surfaces
    visualize =1;
    if visualize
        figure('Name','Data Conversion','NumberTitle','off','color','w');
        title('Surface data');
        patch('Faces',lhfaces+1,'Vertices',lhvtx,'FaceColor','interp','EdgeColor','none','Facevertexcdata',lhix)
        patch('Faces',rhfaces+1,'Vertices',rhvtx,'FaceColor','interp','EdgeColor','none','Facevertexcdata',rhix)
        axis equal
        view([-90,0])
        axis off
        rotate3d on
    end
end

% close all
bnd_L = [];
bnd_L.pnt = lhvtx;
bnd_L.tri = lhfaces+1;
bnd_L.idx = lhix;

bnd_R = [];
bnd_R.pnt = rhvtx;
bnd_R.tri = rhfaces+1;
bnd_R.idx = rhix;


ldat=[lhix;rhix];
if flag.plotatlas ==1
    
    figure;
    ft_plot_mesh(bnd_L, 'faecolor', 'brain',  'vertexcolor', ...
        lhix, 'facealpha', 1);
    hold on
    ft_plot_mesh(bnd_R, 'faecolor', 'brain',  'vertexcolor', ...
        rhix, 'facealpha', 1);
    title('HCP-MMP1, 180 ROIs')
    
    %
    % bnd_LR = [];
    % bnd_LR.pnt = [lhvtx;rhvtx];
    % bnd_LR.tri = [lhfaces+1;rhfaces+1];
    %
    % figure;
    % ft_plot_mesh(bnd_LR, 'faecolor', 'brain',  'vertexcolor', ...
    % [lhix;2*rhix], 'facealpha', 1);
end