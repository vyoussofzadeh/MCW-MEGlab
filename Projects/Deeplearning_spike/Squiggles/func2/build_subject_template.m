function T = build_subject_template(Xcell_spike, fs, varargin)
% BUILD_SUBJECT_TEMPLATE  Build a per-subject spike template from TRAIN spikes only.
%   T = build_subject_template(Xcell_spike, fs, 'trim_ms', 250)
%
% Inputs
%   Xcell_spike : 1xN cell, each [C x T_i] epoch from the SAME subject (Spike only, TRAIN split)
%   fs          : sampling rate (Hz)
% Params
%   trim_ms     : half-window around GFP peak (ms). Output width Tw = 2*trim_ms + 1 (samples)
%
% Output struct T:
%   .wave  : [C x Tw] median waveform template (peak-centered)
%   .topo  : [C x 1]  RMS spatial map (scaled)
%   .cov   : [C x C]  channel covariance of template
%   .gfp   : [1 x Tw] GFP of template
%   .fs    : sampling rate (Hz)

p = inputParser; addParameter(p,'trim_ms',250,@(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(p,varargin{:});
trim = max(1, round(p.Results.trim_ms * fs / 1000));   % samples

assert(iscell(Xcell_spike) && ~isempty(Xcell_spike), 'Xcell_spike must be a non-empty cell array.');

C = size(Xcell_spike{1},1);

% 1) Align each spike to its own GFP peak; center to median peak index
gfps = cellfun(@(X) std(single(X),[],1), Xcell_spike, 'UniformOutput', false);
pidx = cellfun(@(g) find(g==max(g),1,'first'), gfps);           % peak index per epoch
c_ref = round(median(double(pidx)));                             % global center

aligned = cell(size(Xcell_spike));
for i = 1:numel(Xcell_spike)
    X  = single(Xcell_spike{i});                                 % [C x T_i]
    sh = pidx(i) - c_ref;                                        % shift so peaks align
    aligned{i} = circshift(X, [0 -sh]);
end

% 2) Safe trim+pad to fixed peak-centered window Tw
Tw = 2*trim + 1;                                                 % desired width
alignedFixed = cell(size(aligned));
for i = 1:numel(aligned)
    A  = aligned{i}; Ti = size(A,2);
    c0 = 1 + mod(c_ref-1, Ti);                                   % valid center in 1..Ti
    L1i = max(1, c0 - trim);  L2i = min(Ti, c0 + trim);
    S = A(:, L1i:L2i);                                           % [C x w]
    w = size(S,2);
    out = zeros(C, Tw, 'like', A);
    idxCenter = (c0 - L1i) + 1;                                  % center of S (1-based)
    desiredCent = trim + 1;
    startCol = desiredCent - (idxCenter - 1);
    startCol = max(1, min(Tw - w + 1, startCol));
    out(:, startCol:startCol+w-1) = S;
    alignedFixed{i} = out;                                       % [C x Tw]
end
aligned = alignedFixed;

% 3) Template stats
W    = cat(3, aligned{:});                                       % [C x Tw x N]
wave = median(W, 3, 'omitnan');                                  % [C x Tw]
topo = rms(wave, 2); topo = topo ./ (median(abs(topo)) + eps);   % [C x 1]
covT = cov(wave.');                                              % [C x C]
gfpT = std(wave, [], 1);                                         % [1 x Tw]

T = struct('wave',wave,'topo',topo,'cov',covT,'gfp',gfpT,'fs',fs);
end
