function B = feats_spectral_moving(X, fs, bands, winSec, hopSec)
% X: C×T  -> struct with fields bp1..bpK, each C×T
if nargin<3 || isempty(bands)
    bands = [0.5 4; 4 8; 8 13; 13 30; 30 70];  % d ? a ß ?
end
if nargin<4 || isempty(winSec), winSec = 0.25; end
if nargin<5 || isempty(hopSec), hopSec = 0.05; end

[C,T] = size(X);
if T==0, error('feats_spectral_moving: empty epoch'); end

% window/hop (clamped to T)
wSamp = min(max(1, round(winSec*fs)), T);
hSamp = max(1, round(hopSec*fs));

% frames (>=1)
nFrames  = max(1, floor((T - wSamp)/hSamp) + 1);
idxStart = (0:(nFrames-1))*hSamp + 1;
idxEnd   = min(idxStart + wSamp - 1, T);

K = size(bands,1);
BP = zeros(C, nFrames, K, 'single');                 % C×frames×K

% per-frame bandpowers (relative)
for j = 1:nFrames
    xw = X(:, idxStart(j):idxEnd(j));                % C×winLen
    winLen = size(xw,2);

    % simple within-frame PSD
    nfft = max(128, 2^nextpow2(winLen));
    [Pxx,f] = pwelch(xw.', hamming(winLen), 0, nfft, fs, 'psd');  % F×C
    Pxx = Pxx.'; 
    if numel(f) < 2
        % degenerate spectrum -> zeros
        continue;
    end
    df  = f(2)-f(1) + eps; 
    tot = sum(Pxx,2)*df + eps;

    for k = 1:K
        m = f>=bands(k,1) & f<bands(k,2);
        if any(m)
            BP(:,j,k) = (sum(Pxx(:,m),2)*df) ./ tot;             % C×1
        end
    end
end

% Map frames ? sample grid
tFrames = idxStart + (wSamp-1)/2;          % 1×nFrames
tFull   = (1:T).';                          % T×1

for k = 1:K
    % V is nFrames×C (no squeeze!)
    V = permute(BP(:,:,k), [2 1 3]);       % (frames × C)
    if nFrames == 1
        % single frame: broadcast across time
        mapK = repmat(V, T, 1);            % T×C
    else
        % linear interp along frame axis
        mapK = interp1(tFrames(:), V, tFull, 'linear', 'extrap');  % T×C
    end
    B.(sprintf('bp%d',k)) = mapK.';         % C×T
end
end
