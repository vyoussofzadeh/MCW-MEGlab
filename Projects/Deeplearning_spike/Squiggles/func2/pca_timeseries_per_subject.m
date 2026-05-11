function S = pca_timeseries_per_subject(X, Y, groups, r, fs, subjects)
% X       : {N×1} cell, each C×T (already normalized + fixed-length)
% Y       : N×1 categorical (e.g., {'NoSpike','Spike'})
% groups  : {N×1} cell (subject/session IDs)
% r       : # channel PCA comps to keep (e.g., 32)
% fs      : sampling rate in Hz (or [] to plot in samples)
% subjects: (optional) cellstr of subject IDs to restrict; default = all
%
% Returns struct S with one field per subject (S.(id)) containing:
%   .W (C×r), .explained (r×1 %), .counts (table), .fisher12 (scalar)

if nargin<6 || isempty(subjects), subjects = unique(groups,'stable'); end
if nargin<5, fs = []; end
Kplot = 3;                          % plot first 3 PCs
S = struct;

for s = 1:numel(subjects)
    sid = subjects{s};
    idx = strcmp(groups, sid);
    if ~any(idx), continue; end

    Xs = X(idx); Ys = Y(idx);
    C  = size(Xs{1},1); T = size(Xs{1},2);
    r_ = min([r C]);

    % ---- subject PCA from average covariance ----
    covm = zeros(C,'double'); nE = 0;
    for i=1:numel(Xs)
        Xi = double(Xs{i}); covm = covm + (Xi*Xi.')/size(Xi,2); nE=nE+1;
    end
    covm = covm / max(nE,1);
    [U,Sv] = svd(covm,'econ');
    W      = U(:,1:r_);
    expl   = 100*diag(Sv)/sum(diag(Sv));
    expl   = expl(1:r_);

    % ---- project to PCA time-series (r×T) ----
    projOne = @(Z) W.'*double(Z);
    compAll = cellfun(projOne, Xs, 'UniformOutput', false);

    % split by class
    comp0 = compAll(Ys=='NoSpike');
    comp1 = compAll(Ys=='Spike');

    % mean±SEM
    meanSem = @(Ck) deal( ...
        squeeze(mean(cat(3,Ck{:}),3,'omitnan')), ...
        squeeze(std( cat(3,Ck{:}),0,3,'omitnan')) / sqrt(max(1,numel(Ck))) );
    [m0,s0] = meanSem(comp0);   % r×T
    [m1,s1] = meanSem(comp1);   % r×T

    % simple Fisher separation on PC1PC2 averages
    % build per-epoch 2D features by log-variance of first 2 PCs
    lvar = @(Ck) cell2mat(cellfun(@(a) log(var(a(1:2,:),0,2)+eps).', Ck, 'UniformOutput', false));
    F0 = lvar(comp0); F1 = lvar(comp1);
    mu0 = mean(F0,1); mu1 = mean(F1,1);
    S0  = cov(F0);    S1  = cov(F1);
    fisher12 = sum((mu1-mu0).^2) / max(trace(S0+S1),eps);

    % plot
    if isempty(fs), t = 1:T; xlab='Samples';
    else,            t = (0:T-1)/fs; xlab='Time (s)'; end

    figure('Name', sprintf('Subject %s  PCA time-series (first %d PCs)', sid, min(Kplot,r_)), 'Color','w'); clf
    for k=1:min(Kplot,r_)
        subplot(min(Kplot,r_),1,k); hold on
        shaded(t, m0(k,:), s0(k,:), [0.3 0.5 1.0], 0.15);
        shaded(t, m1(k,:), s1(k,:), [1.0 0.5 0.2], 0.15);
        ylabel(sprintf('PC%d',k));
        if k==1
            title(sprintf('Subject %s | chan-PCA(%d)  top3 explained=%.1f/%.1f/%.1f%%', ...
                  sid, r_, expl(1), expl(min(2,end)), expl(min(3,end))));
            legend({'NoSpike','Spike'}, 'Location','best');
        end
        if k==min(Kplot,r_), xlabel(xlab); end
        grid on; xlim([t(1) t(end)]);
    end
    hold off

    % counts table
    cats  = categories(Ys);
    cts   = countcats(Ys);
    countsTbl = table(cats, cts, 'VariableNames', {'Class','Count'});

    % store
    S.(matlab.lang.makeValidName(sid)).W          = W;
    S.(matlab.lang.makeValidName(sid)).explained  = expl(:);
    S.(matlab.lang.makeValidName(sid)).counts     = countsTbl;
    S.(matlab.lang.makeValidName(sid)).fisher12   = fisher12;
end
end

function shaded(x, m, s, col, a)
fill([x, fliplr(x)], [m+s, fliplr(m-s)], col, 'FaceAlpha', a, 'EdgeColor','none');
plot(x, m, 'Color', col, 'LineWidth', 1.4);
end
