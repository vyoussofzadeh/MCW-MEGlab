function runDetect(~,~)
% Gather GUI parameters → cfg
cfg           = struct();
cfg.dtype     = 'meg';                 % or 'eeg' – up to you
cfg.fs        = fs;
cfg.spikeband = sort([ui.spLo.Value ui.spHi.Value]);   % [low high]
cfg.zthresh   = ui.spThr.Value;
cfg.minterval = 0.03;                  % s   (can expose in GUI later)
cfg.statswin  = 2;                     % s
cfg.revspikes = 0;                     % no interactive review
cfg.statswinType = 'global';           % default

% --- read the WHOLE file (all channels, all samples) -------------
mydat = ft_read_data(filename, 'header', hdr);   % M×N
[~, allSpikes] = detMEGspikes(mydat, labels, cfg);

% sample indices → seconds
ui.spikeT = allSpikes / fs;

% update table & repaint
ui.spikeTable.Data = num2cell(ui.spikeT(:));
updatePlot;

msgbox(sprintf('%d spikes detected in entire file', numel(ui.spikeT)), ...
    'Spike detection');
end
