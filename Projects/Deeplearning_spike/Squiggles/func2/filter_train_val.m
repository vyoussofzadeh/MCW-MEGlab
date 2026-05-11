function [Xtr_f, Xva_f] = filter_train_val(Xtr, Xva, fs, bp, bpOrder, notchHz, notchHarm, notchQ)
% FILTER_TRAIN_VAL  Zero-phase band-pass (+ optional notch) for {C×T} cell arrays.

% ---- Validate & defaults ----
assert(iscell(Xtr) && iscell(Xva), 'Xtr/Xva must be cell arrays of C×T matrices.');
if nargin < 8 || isempty(notchQ),    notchQ = 30; end
if nargin < 7 || isempty(notchHarm), notchHarm = []; end
bp = bp(:).'; assert(numel(bp)==2 && bp(1)>0 && bp(2) < fs/2, 'bp must be [low high] Hz and < Nyquist');
bpOrder = max(1, round(bpOrder));

% ---- Band-pass design ----
[bBP, aBP] = butter(bpOrder, bp/(fs/2), 'bandpass');

% ---- Notch design (optional) ----
useNotch = ~isempty(notchHz) && isnumeric(notchHz) && isscalar(notchHz) ...
         && isfinite(notchHz) && notchHz > 0 && notchHz < fs/2;

if useNotch
    freqs = [notchHz, notchHz .* notchHarm(:)'];
    freqs = freqs(freqs > 0 & freqs < fs/2);
    if isempty(freqs)
        useNotch = false;
        bNZ = 1; aNZ = 1;
    else
        [bNZ, aNZ] = notch_bank_design(freqs, notchQ, fs);
    end
else
    bNZ = 1; aNZ = 1;
end

% ---- Allocate outputs (always) ----
Xtr_f = cell(size(Xtr));
Xva_f = cell(size(Xva));

% ---- Filter helpers ----
    function Y = filt_epoch(X)
        if isempty(X), Y = X; return; end
        Xd = double(X);
        % band-pass
        Y  = filtfilt(bBP, aBP, Xd.').';
        % notch (if enabled or identity TF if disabled)
        if useNotch
            Y = filtfilt(bNZ, aNZ, Y.').';
        end
        % cast back to original type
        Y = cast(Y, 'like', X);
    end

% ---- Filter TRAIN ----
for i = 1:numel(Xtr)
    Xi = Xtr{i};
    if ~isempty(Xi) && ~ismatrix(Xi)
        error('Each epoch must be C×T; got size %s at Xtr{%d}', mat2str(size(Xi)), i);
    end
    Xtr_f{i} = filt_epoch(Xi);
end

% ---- Filter VAL ----
for i = 1:numel(Xva)
    Xi = Xva{i};
    if ~isempty(Xi) && ~ismatrix(Xi)
        error('Each epoch must be C×T; got size %s at Xva{%d}', mat2str(size(Xi)), i);
    end
    Xva_f{i} = filt_epoch(Xi);
end
end

function [b, a] = notch_bank_design(freqs, Q, fs)
% Cascade multiple IIR notches into one DF1 filter
b = 1; a = 1;
for f0 = freqs(:).'
    w0 = f0/(fs/2);                        % normalized
    BW = w0 / Q;                           % bandwidth from Q
    [bi, ai] = iirnotch(w0, BW);
    b = conv(b, bi); a = conv(a, ai);
end
end
