function M = feats_morph_time(X, fs)
% X: C×T; returns struct with C×T (maps) and C×1 (scalars)
C = size(X,1); T = size(X,2);

% derivatives
M.d1  = [diff(X,1,2), X(:,end)];
M.d2  = [diff(X,2,2), X(:,end), X(:,end)];

% line length in a short window (20 ms)
wLL = max(1, round(0.02*fs));
M.ll = movsum(abs(diff(X,1,2)), wLL, 2);       % C×(T-1) -> ragged edge OK

% Teager-Kaiser energy
Xp   = [X(:,1), X, X(:,end)];
M.tke = X.^2 - Xp(:,1:end-2).*Xp(:,3:end);     % C×T

% Hilbert envelope power
M.env = abs(hilbert(double(X).').').^2;        % C×T

% zero-crossing rate (per channel scalar)
M.zcr = sum(abs(diff(sign(X),1,2))==2,2) / T;  % C×1
end
