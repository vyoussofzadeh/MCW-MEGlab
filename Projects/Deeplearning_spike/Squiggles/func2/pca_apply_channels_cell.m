function [Zcells] = pca_apply_channels_cell(Xcells, W, mC, K)
% Transform each epoch C×Tfix -> K×Tfix using top-K PCs
Zcells = cell(size(Xcells));
for i=1:numel(Xcells)
    X = single(Xcells{i});
    X0 = X - mC;                   % center with TRAIN mean
    Z  = W(:,1:K)' * X0;           % K×Tfix
    Zcells{i} = Z;
end
end
