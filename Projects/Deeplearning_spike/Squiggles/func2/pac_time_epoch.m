function PAC = pac_time_epoch(X, fs, phaseBands, ampBands, nCycles, hopSec, nSurr)
% Time-resolved PAC per channel for one epoch.
% X: C×T; fs: Hz
% phaseBands: P×2 [low high] Hz (e.g., [1 4; 4 8])
% ampBands  : A×2 [low high] Hz (e.g., [13 30; 30 55])
% nCycles   : #cycles of the *low* band per bin (default 3)
% hopSec    : hop between bins (s), default 0.01 (10 ms)
% nSurr     : #surrogates for zPAC (default 100)

if nargin<5 || isempty(nCycles), nCycles=3; end
if nargin<6 || isempty(hopSec),  hopSec=0.01; end
if nargin<7 || isempty(nSurr),   nSurr=100; end

[C,T] = size(X);
P = size(phaseBands,1);
A = size(ampBands,1);
hopSamp = max(1, round(hopSec*fs));

% Precompute analytic signals for all bands (whole epoch)
Xp = cell(P,1);  % slow bands analytic
for ip=1:P
    Xp{ip} = band_hilbert(X, fs, phaseBands(ip,:)); % C×T complex
end
Xa = cell(A,1);  % fast bands analytic
for ia=1:A
    Xa{ia} = band_hilbert(X, fs, ampBands(ia,:));   % C×T complex
end

% Bin edges per phase band (bin length depends on its low edge)
bins = cell(P,1); centers = cell(P,1);
for ip=1:P
    flo = phaseBands(ip,1);
    winSec = max(nCycles/flo, 0.04);                 % =nCycles, min 40 ms
    wSamp  = min(max(1, round(winSec*fs)), T);       % clamp to T
    nb     = max(1, floor((T - wSamp)/hopSamp) + 1); % =1 bin
    sIdx   = (0:(nb-1))*hopSamp + 1;
    eIdx   = min(sIdx + wSamp - 1, T);
    bins{ip}    = [sIdx(:) eIdx(:)];                 % nb×2
    centers{ip} = round((sIdx + eIdx)/2);            % 1×nb (indices)
end

% Allocate outputs
PAC.tbins = centers;            % cell{ip} of 1×nb indices
PAC.z     = cell(P,A);          % zPAC C×nb
PAC.raw   = cell(P,A);          % raw MVL C×nb
PAC.mu0   = cell(P,A);          % surrogate mean C×nb

% Loop bands and bins
for ip=1:P
    phi_all = angle(Xp{ip});                    % C×T
    binsPA  = bins{ip}; nb = size(binsPA,1);
    for ia=1:A
        Aenv_all = abs(Xa{ia});                 % C×T
        raw  = zeros(C, nb, 'single');
        mu0  = zeros(C, nb, 'single');
        sig0 = zeros(C, nb, 'single');

        % Precompute a few random lags for surrogates
        % (lagging envelope breaks phase?amp relation)
        sLags = randi([0 T-1], nSurr, 1);

        for b = 1:nb
            s = binsPA(b,1); e = binsPA(b,2);
            phi = phi_all(:, s:e);              % C×w
            Aev = Aenv_all(:, s:e);             % C×w

            % raw MVL per channel
            raw(:,b) = abs(mean(Aev .* exp(1j*phi), 2));

            % Surrogates: circularly shift Aev along time
            pacs = zeros(C, nSurr, 'single');
            w = e - s + 1;
            for s0 = 1:nSurr
                lag = mod(sLags(s0), w);
                Ash = Aev(:, [lag+1:w, 1:lag]);     % C×w
                pacs(:,s0) = abs(mean(Ash .* exp(1j*phi), 2));
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
