function [W, mC, expl, K] = pca_fit_channels_cell(Xcells, varKeep)
% Fit PCA across channels using TRAIN epochs (Xcells: {N} of C×Tfix)
% Returns:
%   W    : C×C eigenvectors (columns = PCs)
%   mC   : C×1 channel mean used to center
%   expl : variance explained (%) per PC
%   K    : #PCs to reach varKeep (e.g., 95)

if nargin<2, varKeep = 95; end
C = size(Xcells{1},1);
% Concatenate across time & epochs to a big C×M matrix
Xcat = [];  % careful with memory; you can also stream if huge
for i=1:numel(Xcells), Xcat = [Xcat, single(Xcells{i})]; end  %#ok<AGROW>
mC = mean(Xcat, 2);                   % channel mean
X0 = Xcat - mC;                       % center
% SVD for PCA
[U,S,~] = svd(X0, 'econ');            % U: C×r
s = diag(S);                          % singular values
expl = 100 * (s.^2) / sum(s.^2);      % % variance explained per PC
K = find(cumsum(expl) >= varKeep, 1); 
if isempty(K), K = size(U,2); end
W = U;                                 % eigenvectors (columns)
end
