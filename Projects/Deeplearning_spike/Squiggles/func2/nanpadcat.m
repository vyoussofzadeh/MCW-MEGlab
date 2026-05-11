function M = nanpadcat(vecs)
% vecs: 1xN cell of row vectors -> M: [Tmax x N] with NaN padding
Tmax = max(cellfun(@numel, vecs));
N = numel(vecs);
M = nan(Tmax, N, 'like', single(vecs{1}));
for i = 1:N
    v = vecs{i}(:);
    M(1:numel(v), i) = v;
end
end
