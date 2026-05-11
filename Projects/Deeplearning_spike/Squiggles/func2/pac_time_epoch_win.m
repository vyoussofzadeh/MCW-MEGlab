function PAC = pac_time_epoch_win(X, fs, phaseBands, ampBands, winMs, hopMs, nSurr, onsetIdx)
% Time-resolved PAC (MVL) per channel with overlapping bins.
% X         : C×Tfix epoch
% fs        : Hz
% phaseBands: P×2 [low high] Hz (e.g., [1 4; 4 8])
% ampBands  : A×2 [low high] Hz (e.g., [13 30; 30 55])
% winMs     : bin length (ms), e.g., 100
% hopMs     : hop length (ms), e.g., 20
% nSurr     : #surrogates for z-PAC, e.g., 100
% onsetIdx  : sample index of spike onset (default=round(T/2))
%
% Returns:
%   PAC.z{ip,ia}      : C×nbins zPAC per bin
%   PAC.raw{ip,ia}    : C×nbins raw MVL
%   PAC.mu0{ip,ia}    : C×nbins surrogate mean
%   PAC.tbin_sec{ip}  : 1×nbins time (s) relative to onset (0 at spike)
%   PAC.bins{ip}      : nbins×2 [start end] sample indices per bin

if nargin<8 || isempty(onsetIdx), onsetIdx = round(size(X,2)/2); end
if nargin<7 || isempty(nSurr),    nSurr = 100; end

[C,T] = size(X);
P = size(phaseBands,1); A = size(ampBands,1);
wSamp = min(max(1, round((winMs/1000)*fs)), T);   % clamp to epoch length
hSamp = max(1, round((hopMs/1000)*fs));
nb    = max(1, floor((T - wSamp)/hSamp) + 1);
sIdx  = (0:(nb-1))*hSamp + 1;
eIdx  = min(sIdx + wSamp - 1, T);
tbin  = round((sIdx + eIdx)/2);                   % centers
t0rel = (tbin - onsetIdx) / fs;                   % seconds relative to onset

% Analytic signals for all bands (whole epoch)
Xp = cell(P,1); Xa = cell(A,1);
for ip=1:P, Xp{ip} = band_hilbert(X, fs, phaseBands(ip,:)); end
for ia=1:A, Xa{ia} = band_hilbert(X, fs, ampBands(ia,:));   end

PAC.z = cell(P,A); PAC.raw = cell(P,A); PAC.mu0 = cell(P,A);
PAC.tbin_sec = cell(P,1); PAC.bins = cell(P,1);

for ip=1:P
    phi_all = angle(Xp{ip});             % C×T
    PAC.tbin_sec{ip} = t0rel;
    PAC.bins{ip}     = [sIdx(:) eIdx(:)];
    for ia=1:A
        Aenv_all = abs(Xa{ia});          % C×T
        raw  = zeros(C, nb, 'single'); 
        mu0  = zeros(C, nb, 'single'); 
        sig0 = zeros(C, nb, 'single');
        % precompute surrogate lags (per bin we shift within-bin)
        sLags = randi([0 wSamp-1], nSurr, 1);

        for b = 1:nb
            s = sIdx(b); e = eIdx(b); w = e-s+1;
            phi = phi_all(:, s:e);       % C×w
            Aev = Aenv_all(:, s:e);      % C×w

            raw(:,b) = abs(mean(Aev .* exp(1j*phi), 2));

            pacs = zeros(C, nSurr, 'single');
            for k=1:nSurr
                lag = mod(sLags(k), w);
                Ash = Aev(:, [lag+1:w, 1:lag]);   % circular shift within bin
                pacs(:,k) = abs(mean(Ash .* exp(1j*phi), 2));
            end
            mu0(:,b)  = mean(pacs,2);
            sig0(:,b) = std(pacs,0,2) + 1e-9;
        end
        PAC.raw{ip,ia} = raw;
        PAC.mu0{ip,ia} = mu0;
        PAC.z{ip,ia}   = (raw - mu0) ./ sig0;   % zPAC per bin
    end
end
end

function Xh = band_hilbert(X, fs, band)
[b,a] = butter(4, band/(fs/2), 'bandpass');
Xf = filtfilt(b,a,double(X.')).';
Xh = hilbert(Xf.').';
end
