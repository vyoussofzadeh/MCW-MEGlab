function maps = buildSpatialMaps(posXY, T)
% posXY: [C x 2] in any range (will be normalized internally)
% T    : fixed time length after padding/cropping
% out  : [C x T x 2] planes (X-map, Y-map) broadcast across time
posXY = normalize(posXY, "range", [-1 1]);
C = size(posXY,1);
Xmap = repmat(posXY(:,1), 1, T);
Ymap = repmat(posXY(:,2), 1, T);
maps = cat(3, Xmap, Ymap);
end
