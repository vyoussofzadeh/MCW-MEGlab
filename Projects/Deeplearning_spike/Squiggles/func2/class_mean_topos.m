function [Topo, counts] = class_mean_topos(Xcells, Y, fs, win_ms, alignSpikeSign)
% Xcells : {N×1}, each [C×T] epoch (already filtered & scaled)
% Y      : N×1 categorical with categories {'NoSpike','Spike'}
% fs     : sampling rate (Hz)
% win_ms : window length around each epoch's GFP peak (e.g., 80150)
% alignSpikeSign : true/false  flip Spike vectors to a common polarity

if nargin<4 || isempty(win_ms),      win_ms = 120; end
if nargin<5 || isempty(alignSpikeSign), alignSpikeSign = true; end

W = max(1, round(win_ms*fs/1000));
C = size(Xcells{1},1);
spk = []; nsp = [];

for i = 1:numel(Xcells)
    X = Xcells{i};
    g = std(X,[],1); [~,ip] = max(g);
    i1 = max(1, ip - floor(W/2));
    i2 = min(size(X,2), ip + floor(W/2));
    topo = mean(X(:, i1:i2), 2);         % C×1
    if Y(i)=="Spike",   spk(:,end+1) = topo;  %#ok<AGROW>
    else                nsp(:,end+1) = topo;  %#ok<AGROW>
    end
end

% Optional: align spike polarity so opposite-signed spikes don't cancel
if alignSpikeSign && ~isempty(spk)
    ref = spk(:,1);
    for k = 1:size(spk,2)
        if dot(spk(:,k), ref) < 0, spk(:,k) = -spk(:,k); end
    end
end

Topo.Spike   = mean(spk, 2, 'omitnan');   % C×1
Topo.NoSpike = mean(nsp, 2, 'omitnan');   % C×1
counts       = struct('Spike',size(spk,2), 'NoSpike',size(nsp,2));

% Optional zero-center (nice for visualization)
Topo.Spike   = Topo.Spike   - median(Topo.Spike,'omitnan');
Topo.NoSpike = Topo.NoSpike - median(Topo.NoSpike,'omitnan');
end
