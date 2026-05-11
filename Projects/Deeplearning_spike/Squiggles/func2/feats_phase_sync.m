function P = feats_phase_sync(X, fs, neighK, XY)
% X: C×T; neighK: k for k-NN; XY: C×2 (for neighbors)
% Returns:
%   P.ph   : C×T instantaneous phase
%   P.lcv  : C×T local circular variance (0=coherent,1=diverse)
%   P.envC : C×1 envelope-corr node strength in 200 ms window
if nargin<3, neighK=3; end
H  = hilbert(double(X).').';      % C×T complex analytic
ph = angle(H);
P.ph = ph;

% neighbors by distance
idxK = knnsearch(XY, XY, 'K', min(neighK+1, size(XY,1))); % include self first
idxK = idxK(:,2:end); % drop self

% local circular variance across K neighbors
C = size(X,1); T = size(X,2); lcv = zeros(C,T);
for c=1:C
    neigh = idxK(c,:);
    e = exp(1j*ph(neigh,:));             % K×T
    m = mean(e,1);                        % 1×T
    lcv(c,:) = 1 - abs(m);                % C×T
end
P.lcv = lcv;

% envelope correlation node strength over last 200 ms
w = max(1, round(0.2*fs)); E = abs(H); 
R = corr(E(:, end-w+1:end).');          % C×C
P.envC = sum(abs(R),2) - 1;             % C×1
end
