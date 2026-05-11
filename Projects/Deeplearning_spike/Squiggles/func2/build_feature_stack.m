function Xstack = build_feature_stack(X, fs, sens, useSpatial)
% X: C×T (already filtered + scaled + normalized)
% Returns: C×T×D stack with {raw, d1, d2, env, (optional) smooth, edge}

% Time-domain
d1  = [diff(X,1,2), X(:,end)];
d2  = [diff(X,2,2), X(:,end), X(:,end)];
env = abs(hilbert(double(X).').').^2;

if nargin>=3 && useSpatial
    XY = sens.chanpos(:,1:2); % reorder to match Xs refLbl before calling
    D  = squareform(pdist(XY)); W = exp(-(D.^2)/(2*(median(D(:))/3)^2)); W(1:size(X,1)+1:end)=0;
    L  = diag(sum(W,2)) - W;
    lam = 0.1;
    Xs = (eye(size(X,1)) + lam*L) \ X;
    Xe = X - Xs;
    Xstack = cat(3, X, d1, d2, env, Xs, Xe);
else
    Xstack = cat(3, X, d1, d2, env);
end
end
