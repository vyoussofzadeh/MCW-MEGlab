function PAC = pac_time_epoch_fixedwin(X, fs, phaseBands, ampBands, winSec, hopSec, nSurr)
% Time-resolved PAC using explicit window length (winSec) and hop (hopSec).
% Returns PAC.z{ip,ia}: C×nbins; PAC.tbins{ip}: bin centers (samples)

if nargin<6 || isempty(hopSec), hopSec=0.01; end
if nargin<7 || isempty(nSurr), nSurr=100; end
[C,T] = size(X);
P = size(phaseBands,1); A = size(ampBands,1);
hopSamp = max(1, round(hopSec*fs));
wSamp   = min(max(1, round(winSec*fs)), T);     % clamp to T
nb      = max(1, floor((T - wSamp)/hopSamp) + 1);
sIdx    = (0:(nb-1))*hopSamp + 1;
eIdx    = min(sIdx + wSamp - 1, T);

% full-epoch analytic signals for each band (once)
Xp = cell(P,1); Xa = cell(A,1);
for ip=1:P, Xp{ip} = band_hilbert(X, fs, phaseBands(ip,:)); end
for ia=1:A, Xa{ia} = band_hilbert(X, fs, ampBands(ia,:)); end

PAC.tbins = cell(P,1);
PAC.z     = cell(P,A); PAC.raw = cell(P,A); PAC.mu0 = cell(P,A);
lags = randi([0 wSamp-1], nSurr, 1);   % lags within-window

for ip=1:P
    phi_all = angle(Xp{ip});
    tb = round((sIdx + eIdx)/2);
    PAC.tbins{ip} = tb;

    for ia=1:A
        Aenv_all = abs(Xa{ia});
        raw  = zeros(C, nb, 'single'); mu0 = raw; sig0 = raw;

        for b=1:nb
            s=sIdx(b); e=eIdx(b); w=e-s+1;
            phi = phi_all(:,s:e);
            Aev = Aenv_all(:,s:e);

            raw(:,b) = abs(mean(Aev .* exp(1j*phi), 2));

            pacs = zeros(C, nSurr,'single');
            for s0=1:nSurr
                lag = mod(lags(s0), w);
                Ash = Aev(:, [lag+1:w, 1:lag]);      % circular shift inside window
                pacs(:,s0) = abs(mean(Ash .* exp(1j*phi), 2));
            end
            mu0(:,b)  = mean(pacs,2);
            sig0(:,b) = std(pacs,0,2)+1e-9;
        end
        PAC.raw{ip,ia} = raw;
        PAC.mu0{ip,ia} = mu0;
        PAC.z{ip,ia}   = (raw - mu0) ./ sig0;
    end
end
end

function Xh = band_hilbert(X, fs, band)
[b,a] = butter(4, band/(fs/2), 'bandpass');
Xf = filtfilt(b,a,double(X.')).';
Xh = hilbert(Xf.').';
end
