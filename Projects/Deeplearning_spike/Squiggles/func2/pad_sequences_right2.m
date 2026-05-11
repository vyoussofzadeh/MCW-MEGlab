function [Xpad,maxLen] = pad_sequences_right2(X,pct)
% Right-pad/crop to a percentile length (e.g., pct=95)
if nargin<2, pct=95; end
T = cellfun(@(z) size(z,2), X);
maxLen = max(1, round(prctile(T,pct)));
C = size(X{1},1);
Xpad = cell(size(X));
for i=1:numel(X)
    xi = X{i};
    Ti = size(xi,2);
    if Ti > maxLen
        xi = xi(:,1:maxLen);
    elseif Ti < maxLen
        xi = [xi, zeros(C, maxLen-Ti, 'like', xi)];
    end
    Xpad{i} = xi;
end
end

function Ximg = toImageBatch(X, extraMaps)
% Convert {C×T} cell -> 4-D array [C×T×Cin×N]
% extraMaps (optional): [C×T×K] to concatenate per sample
N = numel(X);
C = size(X{1},1);
T = size(X{1},2);
K = 0;
if nargin>1 && ~isempty(extraMaps), K = size(extraMaps,3); end
Ximg = zeros(C,T,1+K,N,'single');
for i=1:N
    Ximg(:,:,1,i) = single(X{i});
    if K>0, Ximg(:,:,2:end,i) = single(extraMaps); end
end
end

function maps = buildSpatialMaps(posXY, T)
% posXY: [C×2] normalized to [-1,1]; returns [C×T×2]
C = size(posXY,1);
Xmap = repmat(posXY(:,1),1,T);
Ymap = repmat(posXY(:,2),1,T);
maps = cat(3, Xmap, Ymap);   % [C×T×2]
end
