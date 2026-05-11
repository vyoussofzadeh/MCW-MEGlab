function [Xtr_s, Xva_s] = smooth_train_val(Xtr, Xva, Fs, cfg)
% Smooth each {C x T} epoch in Xtr/Xva along time.
% Fs  : sampling rate (Hz)
% cfg : struct with fields depending on cfg.type:
%   cfg.type = 'movmean'    -> cfg.win_ms
%            = 'sgolay'     -> cfg.frame_ms, cfg.order (e.g., 3)
%            = 'gauss'      -> cfg.sigma_ms (FWHM2.355*sigma), cfg.win_ms (optional)
%            = 'median'     -> cfg.win_ms
%            = 'lowpass'    -> cfg.fc (Hz), cfg.order (e.g., 4)
%
% All methods are zero-lag (filtfilt or symmetric kernels).

Xtr_s = cellfun(@(x) smooth_epoch(x, Fs, cfg), Xtr, 'UniformOutput', false);
Xva_s = cellfun(@(x) smooth_epoch(x, Fs, cfg), Xva, 'UniformOutput', false);
end

function Xs = smooth_epoch(X, Fs, cfg)
% X: C x T
Xs = X; if isempty(X), return; end
switch lower(cfg.type)
    case 'movmean'
        w = ms2samp(cfg.win_ms, Fs);
        Xs = movmean(X, w, 2, 'Endpoints','shrink');  % symmetric window
    case 'median'
        w = ms2samp(cfg.win_ms, Fs); if mod(w,2)==0, w=w+1; end
        Xs = medfilt1(X.', w, [], 1, 'truncate').';   % filter along time
    case 'sgolay'
        frame = ms2samp(cfg.frame_ms, Fs); if mod(frame,2)==0, frame=frame+1; end
        ord   = min(cfg.order, frame-1);
        Xs    = sgolayfilt(X.', ord, frame, [], 1).'; % no phase lag
    case 'gauss'
        % build symmetric normalized kernel
        sig = cfg.sigma_ms/1000 * Fs;                 % in samples
        if isfield(cfg,'win_ms') && ~isempty(cfg.win_ms)
            W = ms2samp(cfg.win_ms, Fs);
        else
            W = max(5, 2*ceil(4*sig)+1);              % ~±4s
        end
        if mod(W,2)==0, W=W+1; end
        k = gausswin(W, 1/(2*sig^2));                 % approximates exp(-n^2/(2s^2))
        k = k / sum(k);
        Xs = conv2(X, k.', 'same');                   % 0-lag symmetric conv
    case 'lowpass'
        fc = cfg.fc; ord = getfield(cfg,'order',4); %#ok<GFLD>
        [b,a] = butter(ord, fc/(Fs/2), 'low');
        Xs = filtfilt(b, a, double(X.')).';
    otherwise
        error('Unknown smoothing type: %s', cfg.type);
end
end

function w = ms2samp(win_ms, Fs), w = max(1, round(win_ms/1000 * Fs)); end
