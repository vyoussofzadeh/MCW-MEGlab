function detectSpikes(bp, zThr)
if nargin<1, bp  = [5 50]; end
if nargin<2, zThr = 12;    end

rngBeg = 1;
rngEnd = nsamples;

chIdx = getCurrentIdx();
if isempty(chIdx), return; end

raw = ft_read_data(filename,'header',hdr, ...
    'begsample',rngBeg,'endsample',rngEnd, ...
    'chanindx',chIdx);

% simple peak-finder:
[b,a] = butter(4, bp/(fs/2));
fil   = filtfilt(b,a, raw')';          % band-pass
zsig  = fil./std(fil(:),1,'omitnan');
[rows,cols] = find(zsig > zThr);

% absolute times (s)
spikeTimes = (rngBeg-1 + cols)'/fs;

ui.spikeT = unique(spikeTimes);   % store
refreshSpikeTable();              % <-- add this line
updatePlot();
msgbox(sprintf('%d spikes detected',numel(ui.spikeT)),'Detect spikes');
end
