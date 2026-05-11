function y = aug_mix_frequency(x, partner, aug)
% Mix magnitude spectra of x with a partner (keep x phase); channelwise
T = min(size(x,2), size(partner,2));
if T <= 1, y = x; return; end

x = x(:,1:T); p = partner(:,1:T);
X = fft(x, [], 2); P = fft(p, [], 2);
alpha = aug.mix_alpha_min + (aug.mix_alpha_max-aug.mix_alpha_min)*rand;
mag = (1-alpha)*abs(X) + alpha*abs(P);
phs = angle(X);
Z = mag .* exp(1j*phs);
y = real(ifft(Z, [], 2));
end
