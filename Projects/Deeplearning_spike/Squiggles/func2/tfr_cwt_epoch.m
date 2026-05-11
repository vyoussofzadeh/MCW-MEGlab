function [Pow, f] = tfr_cwt_epoch(X, fs, fminmax)
% X: C×T, returns Pow: C×F×T (power), f: F×1 frequencies
if nargin<3, fminmax = [2 90]; end
C = size(X,1);
% Choose a compact frequency grid (log-spaced)
f = logspace(log10(fminmax(1)), log10(fminmax(2)), 40);  % 40 freqs ? good balance
Pow = zeros(C, numel(f), size(X,2), 'single');
for c = 1:C
    % cwtfilterbank: good time-freq tradeoff
    fb = cwtfilterbank('SignalLength',size(X,2),'SamplingFrequency',fs, ...
                       'FrequencyLimits',fminmax,'VoicesPerOctave',12, ...
                       'TimeBandwidth',60);
    % Get full spectrum; interpolate to f grid (optional)
    [wt, freqs] = cwt(double(X(c,:)), 'FilterBank', fb);   % wt: F0×T
    % Interp to your f grid so all channels share the same F axis
    Pow(c,:,:) = single(abs(interp1(freqs, wt, f, 'linear','extrap')).^2);
end
end
