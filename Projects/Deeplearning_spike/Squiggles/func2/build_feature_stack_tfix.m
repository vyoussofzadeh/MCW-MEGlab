function Xstack = build_feature_stack_tfix(X, fs, sens, useSpatial)
% Build C×T×D stack for one epoch X (C×Tfix).
% D = 4: {raw, d1, d2, envelope};  D = 6 if useSpatial: + {spatial-smooth, spatial-edge}
% - X must already be C×Tfix (fixed length)
% - If useSpatial=true, pass FieldTrip 'sens' aligned to X's channel order (refLbl)

C = size(X,1); T = size(X,2);

% ---- time-domain features ----
d1  = [diff(X,1,2), X(:,end)];                    % onset slope
d2  = [diff(X,2,2), X(:,end), X(:,end)];          % curvature
env = abs(hilbert(double(X).').').^2;             % envelope power (C×T)

if nargin>=3 && ~isempty(sens) && useSpatial
    % Build a quick Laplacian from chanpos; sens must be aligned to X's rows first
    XY = sens.chanpos(:,1:2);
    Dm = squareform(pdist(XY));
    sig = median(Dm(:))/3;
    W  = exp(-(Dm.^2)/(2*(sig^2))); 
    W(1:C+1:end) = 0;
    L  = diag(sum(W,2)) - W;                       % graph Laplacian
    lam = 0.1;                                     % smoothness knob
    Xs = (eye(C) + lam*L) \ X;                     % spatially smoothed
    Xe = X - Xs;                                   % spatial high-pass residual
    Xstack = single(cat(3, X, d1, d2, env, Xs, Xe));    % C×T×6
else
    Xstack = single(cat(3, X, d1, d2, env));            % C×T×4
end
end
