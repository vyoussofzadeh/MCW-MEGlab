function ft_topo_mag_window(X, fs, refLbl, win_ms)
% X: [C x T] epoch, refLbl: cellstr MEG labels in row order
% win_ms: [t1 t2] in ms relative to epoch start (e.g., [80 120])

% 1) Build a minimal timelock struct (single average vector over window)
t = (0:size(X,2)-1)/fs;
t1 = win_ms(1)/1000; t2 = win_ms(2)/1000;
idx = t>=t1 & t<=t2;
avg = mean(X(:,idx), 2);           % C×1 topography for the window

tl               = [];
tl.label         = refLbl(:);
tl.dimord        = 'chan_time';
tl.time          = 0;              % a single time so ft_topoplotER is happy
tl.avg           = avg(:)';        % 1×C, but ft expects chan×time -> put as row

% 2) Plot
cfg = [];
cfg.parameter = 'avg';
cfg.xlim      = [0 0];             % the only time we have
cfg.layout    = 'neuromag306mag.lay';   % magnetometer layout
cfg.marker    = 'off';
cfg.comment   = 'no';
cfg.colorbar  = 'yes';
cfg.zlim      = 'maxabs';          % symmetric scale
ft_topoplotER(cfg, tl);
title(sprintf('Mag topomap %d%d ms', win_ms(1), win_ms(2)));
end
