function [Xs, Xe] = feats_spatial_graph(X, XY, lam)
% X: C×T; XY: C×2 sensor coords aligned to Xs channel order (refLbl)
% lam: smoothness weight (e.g., 0.1)
if nargin<3, lam = 0.1; end
C = size(X,1);
D  = squareform(pdist(XY));
sigma = median(D(:))/3;
W  = exp(-(D.^2)/(2*sigma^2)); W(1:C+1:end)=0;
L  = diag(sum(W,2)) - W;         % combinatorial Laplacian

Xs = (eye(C) + lam*L) \ X;       % spatially smoothed
Xe = X - Xs;                     % spatial high-pass residual
end

function Xe_abs = feats_spatial_edge_abs(Xe)
Xe_abs = abs(Xe);  % nonlinearity; or use signed edge (Xe) directly
end
