function pca_plot3_spike(Xcell, Y, r)
% 3D PCA plot of epochs using channel-PCA -> log-variance features.
% Xcell : {N×1} cell, each C×T (already normalized & fixed-length)
% Y     : N×1 categorical ('NoSpike','Spike', ...)
% r     : # channel-PCA comps to keep for features (e.g., 32)

assert(~isempty(Xcell) && numel(Xcell)==numel(Y), 'Size mismatch X vs Y');
C = size(Xcell{1},1);

% ----- 1) Channel PCA from average covariance (fit on provided set) -----
covm = zeros(C,'double'); nE = 0;
for i=1:numel(Xcell)
    Xi = double(Xcell{i});                      % C x T
    covm = covm + (Xi*Xi.')/size(Xi,2);
    nE = nE + 1;
end
covm = covm / max(nE,1);
[W,S] = svd(covm,'econ');                       % W: C×C
expl_chan = 100*diag(S)/sum(diag(S));
r = min([r, C]); W = W(:,1:r);                  % take top-r channel PCs

% ----- 2) Epoch-level features: log-variance of channel-PC timecourses -----
feat = @(Xc) cell2mat(cellfun(@(x) log(var((W.'*double(x)),0,2)+eps).', ...
                               Xc, 'UniformOutput', false));     % N x r
F = feat(Xcell);

% Standardize features across epochs (removes scale bias)
muF = mean(F,1); sdF = std(F,0,1) + eps;
Z   = (F - muF) ./ sdF;

% ----- 3) PCA across epochs to 3D -----
[coeff3, score3, expl3] = pca(Z, 'NumComponents', 3);  %#ok<ASGLU>

% ----- 4) 3D scatter (Spike vs NoSpike) -----
idx0 = Y=='NoSpike';
idx1 = Y=='Spike';

figure('Name', sprintf('3D PCA of epochs (chan-PCA=%d -> log-var -> PCA(3))', r));
clf; hold on
s0 = scatter3(score3(idx0,1), score3(idx0,2), score3(idx0,3), 12, [0.25 0.45 1.0], 'filled');
s1 = scatter3(score3(idx1,1), score3(idx1,2), score3(idx1,3), 12, [1.0 0.45 0.15], 'filled');
grid on; axis vis3d
xlabel(sprintf('PC1 (%.1f%%)', expl3(1))); 
ylabel(sprintf('PC2 (%.1f%%)', expl3(2)));
zlabel(sprintf('PC3 (%.1f%%)', expl3(3)));
title(sprintf('channel-PCA (%d) -> log-var -> PCA(3D)', r));
legend([s0 s1], {'NoSpike','Spike'}, 'Location','best');
view(40,22);
hold off
end
