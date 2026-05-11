function Xr = pca_project_apply(X, W)
% Apply previously learned projection W (C×r) to a new set X.
Xr = cellfun(@(z) single(W.'*double(z)), X, 'UniformOutput', false);  % r×T
end