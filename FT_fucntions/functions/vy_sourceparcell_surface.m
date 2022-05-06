function [pow_parcel,rois] = vy_sourceparcell_surface(atlas,ImageGridAmp, sMRI)

rois = [];
pow_parcel = [];
for i=1:length(atlas.Scouts)
    pow_parcel(i) = mean(ImageGridAmp(atlas.Scouts(i).Vertices));
    rois{i} = atlas.Scouts(i).Label;
end