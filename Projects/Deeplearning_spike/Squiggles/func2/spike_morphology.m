function F = spike_morphology(w, fs)
%SPIKE_MORPHOLOGY Compute standard morphology features from a 1xT spike template.
% w : 1xT or Tx1 waveform (will be centered on its main positive peak)
% fs: sampling rate (Hz)

w = double(w(:)');                  % row
T = numel(w);
assert(T>=3,'waveform too short');

% make main peak positive
[~,p] = max(abs(w));
if w(p) < 0, w = -w; end

% basic measures
F.peakAmp = w(p);                   % amplitude at main peak
F.energy  = sum(w.^2)/fs;           % signal energy (a.u.*s)

% half-width (at 50% of peak)
half = 0.5*F.peakAmp;
iL = find(w(1:p) <= half, 1, 'last'); if isempty(iL), iL = p; end
iR = p-1 + find(w(p:end) <= half, 1, 'first'); if isempty(iR), iR = p; end
F.halfWidth_s  = (iR - iL)/fs;
F.halfWidth_ms = 1000*F.halfWidth_s;

% slopes (±10 ms window, clipped to edges)
win = max(1, round(0.01*fs));
i1 = max(1, p - win);
i2 = min(T, p + win);
F.upSlope   = (w(p) - w(i1)) / max(eps, (p - i1)/fs);
F.downSlope = (w(i2) - w(p)) / max(eps, (i2 - p)/fs);   % typically negative

% sharpness (discrete second derivative at peak)
if p>1 && p<T
    F.sharpness = w(p-1) - 2*w(p) + w(p+1);
else
    F.sharpness = NaN;
end

% asymmetry (pre vs post half-amp duration)
F.asymmetry = (p - iL) / max(1, (iR - p));

% decay time to 10% of peak
th = 0.1*F.peakAmp;
jR = p-1 + find(w(p:end) <= th, 1, 'first'); if isempty(jR), jR = T; end
F.decay_ms = 1000*(jR - p)/fs;
end
