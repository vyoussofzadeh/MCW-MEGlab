function S = feats_spectral_epoch(X, fs)
% X: C×T  -> spectral scalars per channel (C×1)
% Returns:
%   bp_* : relative band power for d/?/a/ß/? (C×1)
%   sent : spectral entropy (C×1)

% Welch PSD (250 ms window; 50% overlap)
win = max(64, 2^floor(log2(0.25*fs)));   
nover = round(0.5*win);
nfft = max(256, 2^nextpow2(win));
[Pxx,f] = pwelch(X.', hamming(win), nover, nfft, fs,'psd');   % F×C
Pxx = Pxx.'; 
df = f(2)-f(1)+eps; 
tot = sum(Pxx,2)*df + eps;

rel = @(a,b) sum(Pxx(:, f>=a & f<b), 2)*df ./ tot;  % C×1
S.bp_delta = rel(0.5,4);
S.bp_theta = rel(4,8);
S.bp_alpha = rel(8,13);
S.bp_beta  = rel(13,30);
S.bp_gamma = rel(30,70);

Pg = Pxx ./ tot; Pg(Pg<=0) = eps;
S.sent = -sum(Pg.*log(Pg),2) / log(size(Pxx,2));    % spectral entropy
end
