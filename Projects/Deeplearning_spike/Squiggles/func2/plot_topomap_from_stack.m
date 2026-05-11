function plot_topomap_from_stack(Xstack, fs, sens, refLbl, feat, timewin_sec)
% Xstack: C×T×D (from build_feature_stack)
% fs: sampling rate (Hz)
% sens: FieldTrip sensor struct (e.g., data.grad or a sens file) with fields .label, .chanpos/grad
% refLbl: C×1 cellstr in the SAME order as Xstack
% feat: name or index, e.g. 'raw'|'d1'|'d2'|'env' or D-index
% timewin_sec: [t0 t1] in seconds, e.g. [0.05 0.15]

assert(ndims(Xstack)==3, 'Xstack must be C×T×D');
C = size(Xstack,1); T = size(Xstack,2);
t = (0:T-1)/fs;

% --- pick feature ---
featnames = {'raw','d1','d2','env','Xs','Xe'}; % adjust if you added spatial
if ischar(feat) || isstring(feat)
    fi = find(strcmp(featnames, string(feat)),1);
    assert(~isempty(fi), 'Unknown feat "%s"', feat);
else
    fi = feat;
end
Xfeat = squeeze(Xstack(:,:,fi));   % C×T

% --- reorder sens to match refLbl ---
% (FieldTrip expects data.label and sens.label aligned)
[ok, ia] = ismember(refLbl, sens.label);
assert(all(ok), 'Some refLbl missing in sens.label');
sens_reo = sens;
sens_reo.label = sens.label(ia);
if isfield(sens,'chanpos'), sens_reo.chanpos = sens.chanpos(ia,:); end
if isfield(sens,'coilpos'), sens_reo.coilpos = sens.coilpos(ia,:); end
if isfield(sens,'tra') && ~isempty(sens.tra)
    sens_reo.tra = sens.tra(ia, :);
end

% --- build FT raw-like single-trial struct ---
data = [];
data.label    = refLbl(:);
data.fsample  = fs;
data.time     = {t};
data.trial    = {Xfeat};
% attach sensor geometry (grad/sens)
if isfield(sens_reo,'type') || isfield(sens_reo,'balance')
    data.grad = sens_reo;          % for MEG
else
    data.elec = sens_reo;          % for EEG
end

% --- timelock average (just to use ft_topoplotER interface) ---
cfg = [];
timelock = ft_timelockanalysis(cfg, data);  % produces .avg (C×T), .time, and carries grad

% --- layout for plotting ---
lay = ft_prepare_layout([], timelock);      % auto from grad/elec
cfg = [];
cfg.layout = lay;
cfg.marker = 'on';
cfg.comment = 'xlim';
cfg.zlim = 'maxabs';                        % symmetric color range

% --- choose time window & plot ---
if nargin < 6 || isempty(timewin_sec)
    timewin_sec = [t(round(T*0.4)) t(round(T*0.6))]; % default mid window
end
cfg.xlim = timewin_sec;                     % average over this window
figure('Color','w'); 
ft_topoplotER(cfg, timelock);
title(sprintf('Topomap: %s  (%.0f%.0f ms)', featnames{fi}, 1000*timewin_sec(1), 1000*timewin_sec(2)));
end