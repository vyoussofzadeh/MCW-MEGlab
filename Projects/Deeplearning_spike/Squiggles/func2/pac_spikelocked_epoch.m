function [zPAC, PAC, PAC0] = pac_spikelocked_epoch(X, fs, phaseBand, ampBand, winSamp, nSurr)
% X: C×T epoch (filtered/normalized/padded), spike near the center
% phaseBand=[1 4], ampBand=[30 55], winSamp = samples around spike (e.g., ±75ms at 500Hz -> 75)
% Returns per-channel:
%   PAC  : raw MVL (mean vector length)
%   PAC0 : surrogate mean (null)
%   zPAC : (PAC - mean(surr)) / std(surr)

if nargin<6, nSurr=200; end
C = size(X,1); T = size(X,2); ctr = round(T/2);
i1 = max(1, ctr-winSamp); i2 = min(T, ctr+winSamp);
xw = X(:, i1:i2);                         % C×W

% bandpass -> Hilbert
xp = band_hilbert(xw, fs, phaseBand);    % analytic signal, C×W
xa = band_hilbert(xw, fs, ampBand);      % analytic signal, C×W

phi = angle(xp);                          % slow phase
A   = abs(xa);                            % fast amplitude envelope

% MVL PAC per channel
PAC = abs(mean(A .* exp(1j*phi), 2));     % C×1

% Surrogates: circularly shift A by random lags (break phase?amp relation)
PACsurr = zeros(C, nSurr);
W = size(A,2);
for s = 1:nSurr
    lag = randi([0 W-1],1,1);
    Ash = A(:, mod((0:W-1)+lag, W)+1);
    PACsurr(:,s) = abs(mean(Ash .* exp(1j*phi), 2));
end
mu0 = mean(PACsurr,2); sig0 = std(PACsurr,0,2) + 1e-9;
zPAC = (PAC - mu0) ./ sig0;
PAC0 = mu0;
end

function Xh = band_hilbert(X, fs, band)
% zero-phase IIR -> analytic signal
[b,a] = butter(4, band/(fs/2), 'bandpass');
Xf = filtfilt(b,a,double(X.')).';
Xh = hilbert(Xf.').';
end
