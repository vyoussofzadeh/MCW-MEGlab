function [pow_parcel,rois, pow_parcel_count] = do_sourceparcell_surface(cfg_main)

ImageGridAmp = cfg_main.d_in;
atlas = cfg_main.atlas;
% thre = cfg_main.thre;

% raw value
rois = []; pow_parcel = [];
for i=1:length(atlas.Scouts)
    pow_parcel(i) = mean(ImageGridAmp(atlas.Scouts(i).Vertices));
    rois{i} = atlas.Scouts(i).Label;
end

% Counting active vericies approch
pow_parcel_count = [];
for i=1:length(atlas.Scouts)
    parcelval = (ImageGridAmp(atlas.Scouts(i).Vertices));
    pow_parcel_count{i} = parcelval;
%     length(find(parcelval > thre.*max(parcelval)));
end