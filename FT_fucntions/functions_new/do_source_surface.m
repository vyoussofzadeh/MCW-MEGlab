function [p_val, p_max, rois] = do_source_surface(atlas,ImageGridAmp, ~)

rois = []; p_val = [];
for i=1:length(atlas.Scouts)
    p_val{i} = ImageGridAmp(atlas.Scouts(i).Vertices);
    p_max{i} = max(ImageGridAmp(atlas.Scouts(i).Vertices));
    rois{i} = atlas.Scouts(i).Label;
end