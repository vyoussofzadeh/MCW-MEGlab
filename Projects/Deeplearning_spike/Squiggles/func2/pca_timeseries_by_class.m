function [W, compTS] = pca_timeseries_by_class(Xcell, Y, r, fs)
% Plot first 3 PCA component time-series (mean ± SEM) for Spike vs NoSpike.
% Xcell : {N×1} cell, each C×T (already normalized & fixed-length)
% Y     : N×1 categorical with {'NoSpike','Spike'}
% r     : # channel PCA comps to keep (e.g., 32)  we plot the first 3
% fs    : sampling rate in Hz (optional). If omitted, x-axis is samples.
%
% Returns:
%   W       : C×r PCA projection (channels -> comps)
%   compTS  : struct with fields .Spike and .NoSpike, each {N×1} cell of r×T

if nargin<4 || isempty(fs), fs = []; end
assert(~isempty(Xcell) && numel(Xcell)==numel(Y), 'Size mismatch X vs Y');
C = size(Xcell{1},1); T = size(Xcell{1},2);
r = min([r C]);

% ---- 1) Channel PCA from average covariance (fit on provided set) ----
covm = zeros(C,'double'); nE = 0;
for i = 1:numel(Xcell)
    Xi = double(Xcell{i});
    covm = covm + (Xi*Xi.')/size(Xi,2);
    nE = nE + 1;
end
covm = covm / max(nE,1);
[U,S] = svd(covm,'econ'); %#ok<ASGLU>
W = U(:,1:r);   % C×r

% ---- 2) Project each epoch to PCA space: r×T time-series ----
projOne = @(X) W.'*double(X);   % r×T
compAll = cellfun(projOne, Xcell, 'UniformOutput', false);

% split by class
idx0 = (Y=='NoSpike'); idx1 = (Y=='Spike');
compTS.NoSpike = compAll(idx0);
compTS.Spike   = compAll(idx1);

% ---- 3) Compute mean ± SEM across epochs for each class, first 3 PCs ----
Kplot = min(3, r);
meanSem = @(Ck) deal( ...
    squeeze(mean(cat(3, Ck{:}), 3, 'omitnan')), ...   % r×T mean
    squeeze(std( cat(3, Ck{:}), 0, 3, 'omitnan')) / sqrt(max(1,numel(Ck))) ); % r×T SEM

[m0, s0] = meanSem(compTS.NoSpike);  % r×T
[m1, s1] = meanSem(compTS.Spike);    % r×T

% time vector
if isempty(fs), t = 1:T; xlab = 'Samples';
else, t = (0:T-1)/fs; xlab = 'Time (s)'; end

% ---- 4) Plot ----
figure('Name', sprintf('PCA time-series (first %d comps)', Kplot), 'Color','w'); clf
for k = 1:Kplot
    subplot(Kplot,1,k); hold on
    plot_shaded(t, m0(k,:), s0(k,:), [0.3 0.5 1.0], 0.15);
    plot_shaded(t, m1(k,:), s1(k,:), [1.0 0.5 0.2], 0.15);
    ylabel(sprintf('PC%d', k));
    if k==1, title(sprintf('Channel-PCA (%d comps) ? component time-series (mean±SEM)', r)); end
    if k==Kplot, xlabel(xlab); end
    grid on; xlim([t(1) t(end)]);
    if k==1, legend({'NoSpike','Spike'},'Location','best'); end
end
hold off
end

function plot_shaded(x, m, s, col, alphaVal)
% mean ± SEM shaded line (no toolboxes)
fill([x, fliplr(x)], [m+s, fliplr(m-s)], col, 'FaceAlpha', alphaVal, 'EdgeColor','none');
plot(x, m, 'Color', col, 'LineWidth', 1.5);
end
