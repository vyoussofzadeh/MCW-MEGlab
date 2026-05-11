function X7 = build_longview_features(S, fs)
% S: [C x T] (single/double). fs: sampling rate (Hz).
% Returns X7: [C x T x 7] = cat(raw, topo, amp, meanAmp, slope, halfSlope, sharpness)

C = size(S,1); T = size(S,2);
[topo, amp, meanAmp, slope, halfSlope, sharpness] = deal(zeros(C,T,'like',S));

for c = 1:C
    x = S(c,:);

    % crude extrema (you can refine later)
    [pks, locsMax] = findpeaks(x);
    [mins, locsMin] = findpeaks(-x); mins = -mins;

    L = min(numel(locsMax), numel(locsMin)-1);
    for k = 1:L
        tl = locsMin(k); to = locsMax(k); tr = locsMin(k+1);
        if ~(tl<to && to<tr), continue; end

        % half-height crossings (simple proxies)
        hl = (x(to)+x(tl))/2; hr = (x(to)+x(tr))/2;
        thl = tl + find(x(tl:to) >= hl, 1, 'first') - 1;
        thr = to + find(x(to:tr) >= hr, 1, 'last')  - 1;
        thl = max(tl, min(thl, to));  % clamp
        thr = max(to, min(thr, tr));

        al = x(to) - x(tl);                 % left amplitude
        ar = x(tr) - x(to);                 % right amplitude
        spl  = al / max(1, to - tl);        % left slope
        spr  = ar / max(1, tr - to);        % right slope
        hspl = al / max(1, 2*(to - thl));   % half-slope left
        hspr = ar / max(1, 2*(thr - to));   % half-slope right
        ma   = (al*(tr-to) + ar*(to-tl)) / max(tr - tl, 1);  % mean amp over wave
        ampV = max(abs([al ar]));           % peak-like amplitude

        amp(c, thl:thr)        = ampV;
        meanAmp(c, thl:thr)    = ma;
        slope(c, thl:to)       = spl;   slope(c, to:thr)      = spr;
        halfSlope(c, thl:to)   = hspl;  halfSlope(c, to:thr)  = hspr;

        % topo marker (very simple): +1 at maxima, -1 at adjacent minima
        topo(c, to)   = 1;
        topo(c, [tl tr]) = -1;
    end
end

% --- Long-view normalization (~6 s window) on each feature map ---
winW = max(1, round(6*fs));
normalize = @(F) local_longview_norm(F, winW);

topo       = normalize(topo);
amp        = normalize(amp);
meanAmp    = normalize(meanAmp);
slope      = normalize(slope);
halfSlope  = normalize(halfSlope);
sharpness  = normalize(sharpness);  % currently zeros unless you add a sharper metric

% Stack raw + 6 features
X7 = cat(3, S, topo, amp, meanAmp, slope, halfSlope, sharpness);  % [C x T x 7]
end

function Fz = local_longview_norm(F, winW)
% Row-wise (channel-wise) moving z-score with window winW
mu = movmean(F, winW, 2);                 % [C x T]
sd = movstd( F, winW, 0, 2);              % [C x T]
sd(sd < 1e-6) = 1;
Fz = (F - mu) ./ sd;
end
