function [Xr, W] = pca_project_train(X, r)
% Low-rank spatial projection: learn W (C×r) on TRAIN, apply to TRAIN.
C = size(X{1},1); covm = zeros(C); n=0;
for i=1:numel(X)
    xi = double(X{i});
    covm = covm + (xi*xi.')/size(xi,2); n=n+1;
end
covm = covm/n;
[U,~,~] = svd(covm,'econ');
W = U(:,1:r);                                 % C×r
Xr = cellfun(@(z) single(W.'*double(z)), X, 'UniformOutput', false);  % r×T
end

function Xr = pca_project_apply(X, W)
Xr = cellfun(@(z) single(W.'*double(z)), X, 'UniformOutput', false);
end
