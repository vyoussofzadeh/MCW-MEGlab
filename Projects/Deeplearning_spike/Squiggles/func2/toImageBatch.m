function Ximg = toImageBatch(X, extraMaps)
% Convert cell array of {C x T} epochs into 4D array [C x T x Cin x N].
% X         : 1xN cell, each CxT (single/double)
% extraMaps : (optional) CxT x K planes (e.g., X/Y maps), same C,T for all
% Returns   : Ximg of size [C T (1+K) N], single precision

N = numel(X);
assert(N>=1, 'toImageBatch:Empty', 'X is empty.');

C = size(X{1},1);
T = size(X{1},2);
K = 0;

if nargin>1 && ~isempty(extraMaps)
    assert(ndims(extraMaps)==3, 'extraMaps must be CxTxK.');
    assert(size(extraMaps,1)==C && size(extraMaps,2)==T, ...
        'extraMaps size must match CxT of X.');
    K = size(extraMaps,3);
end

Ximg = zeros(C, T, 1+K, N, 'single');

for i = 1:N
    Xi = X{i};
    assert(ismatrix(Xi) && size(Xi,1)==C && size(Xi,2)==T, ...
        'All X{i} must be CxT with same sizes. Consider padding first.');
    Ximg(:,:,1,i) = single(Xi);
    if K>0
        Ximg(:,:,2:end,i) = single(extraMaps);
    end
end
end
