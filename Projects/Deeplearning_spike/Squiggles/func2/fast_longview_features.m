function X = fast_longview_features(S, fs)
% S: [C x T] single/double. Returns X: [C x T x 5]
% Features: raw, moving peak-to-peak, moving mean |x|, moving |dx|, moving |d2x|
S = single(S);
[C,T] = size(S);

% Tunable window sizes (keep small)
wAmp   = max(1, round(0.050*fs));  % ~50 ms
wMean  = max(1, round(0.040*fs));  % ~40 ms
wSlope = max(1, round(0.030*fs));  % ~30 ms

% Differences (pad to T)
d1 = diff(S,1,2);              d1 = padarray(d1,[0 1],'replicate','post');
d2 = diff(S,2,2);              d2 = padarray(d2,[0 2],'replicate','post');

% Vectorized moving ops along time (dim 2)
amp_pp   = movmax(S, wAmp, 2) - movmin(S, wAmp, 2);
meanAbs  = movmean(abs(S), wMean, 2);
slopeAbs = movmean(abs(d1), wSlope, 2);
sharp2   = movmean(abs(d2), wSlope, 2);

% Long-view (6 s) local z-score per channel
winLV = max(1, round(6*fs));
normlv = @(F) (F - movmean(F,winLV,2)) ./ max(movstd(F,winLV,0,2), 1e-6);
amp_pp   = normlv(amp_pp);
meanAbs  = normlv(meanAbs);
slopeAbs = normlv(slopeAbs);
sharp2   = normlv(sharp2);

X = cat(3, S, amp_pp, meanAbs, slopeAbs, sharp2);  % [C x T x 5]
end
