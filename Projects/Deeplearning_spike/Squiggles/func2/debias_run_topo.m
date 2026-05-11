function Xdb = debias_run_topo(X, runMeanTopo)
% Subtract run-mean topomap (C×1) from each time slice
Xdb = X - runMeanTopo;
end

function snrMap = local_snr_map(X, baseIdx)
% X: C×T; baseIdx: indices for a baseline segment
mu = mean(abs(X), 2);
sig = std(X(:, baseIdx), 0, 2) + 1e-6;
snrMap = (mu ./ sig) * ones(1, size(X,2));   % broadcast C×1 to C×T
end
