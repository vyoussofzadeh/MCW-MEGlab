function f = feat_vs_template(X, T)
% FEAT_VS_TEMPLATE  Compare epoch X to subject template T; return 1x6 features.
% X: [C x T] epoch (same subject). T: struct from build_subject_template.
% Output f: [1 x 6] = [r_mean, r_max, r_med, r_topo, r_gfp, d_cov]

X = single(X); C = size(X,1);

% 1) Align by epoch's own GFP peak to template center (no xcorr)
gX = std(X,[],1);
[~, iPeak] = max(gX);
desiredCenter = round((numel(T.gfp)+1)/2);   % trim+1
sh = iPeak - desiredCenter;
X = circshift(X, [0 -sh]);

% 2) Crop/pad to template width Tw with peak centered
Tw = size(T.wave,2); Ti = size(X,2);
Ttrim = floor((Tw-1)/2);
c0 = min(max(desiredCenter,1), Ti);
L1 = max(1, c0 - Ttrim); L2 = min(Ti, c0 + Ttrim);
S  = X(:, L1:L2);
w = size(S,2); out = zeros(C, Tw, 'like', X);
idxCenter = (c0 - L1) + 1;
startCol  = (Ttrim+1) - (idxCenter - 1);
startCol  = max(1, min(Tw - w + 1, startCol));
out(:, startCol:startCol+w-1) = S;
Xw = out;                                     % [C x Tw]

% 3) Similarities
% (a) waveform corr per channel -> aggregates
num = sum(T.wave .* Xw, 2);
den = sqrt(sum(T.wave.^2,2).*sum(Xw.^2,2) + eps);
r_ch   = num ./ (den + eps);
r_mean = mean(r_ch, 'omitnan'); r_max = max(r_ch); r_med = median(r_ch, 'omitnan');

% (b) spatial/topographic correlation (RMS across time)
topoX = rms(Xw, 2); topoX = topoX / (median(abs(topoX)) + eps);
r_topo = corr(topoX, T.topo, 'type','Spearman', 'rows','pairwise');

% (c) covariance distance (log-Euclidean)
covX = cov(Xw.');
[V1,D1] = eig(T.cov + 1e-6*eye(C)); L1 = V1*log(max(D1,1e-9))*V1';
[V2,D2] = eig(covX + 1e-6*eye(C));  L2 = V2*log(max(D2,1e-9))*V2';
d_cov = norm(L1 - L2, 'fro');

% (d) GFP correlation (now equal length)
r_gfp = corr(std(Xw,[],1)', T.gfp', 'rows','pairwise');

% f = [r_mean, r_max, r_med, r_topo, r_gfp, d_cov];


% ... your existing body ...
f = [r_mean, r_max, r_med, r_topo, r_gfp, d_cov];
f = reshape(double(f), 1, []);   % <-- ALWAYS 1×6 row

end
