function [Topo, counts] = group_topo_mean(Xcells, Y, fs, refLbl, win_ms, alignSpikes)
% Xcells : {N×1}  each C×T epoch (already filtered / scaled)
% Y      : N×1 categorical with categories like {'NoSpike','Spike'}
% fs     : sampling rate (Hz)
% refLbl : C×1 channel labels in row order
% win_ms : window length (e.g., 150)
% alignSpikes : true/false  flip spike maps to align polarity

if nargin<5 || isempty(win_ms), win_ms = 150; end
if nargin<6, alignSpikes = true; end
W = max(1, round(win_ms*fs/1000));

C = size(Xcells{1},1);
S = []; N = [];

% 1) collect per-epoch topomaps (mean over window around GFP peak)
for i = 1:numel(Xcells)
    X = Xcells{i};
    g = std(X,[],1); [~, ip] = max(g);
    i1 = max(1, ip - floor(W/2));
    i2 = min(size(X,2), ip + floor(W/2));
    topo = mean(X(:, i1:i2), 2);          % C×1
    if Y(i)=="Spike"
        S(:,end+1) = topo;                %#ok<AGROW>
    else
        N(:,end+1) = topo;                %#ok<AGROW>
    end
end

% 2) optional: align spike polarity so dipoles dont cancel
if alignSpikes && ~isempty(S)
    ref = S(:,1);
    for k = 1:size(S,2)
        if dot(S(:,k), ref) < 0
            S(:,k) = -S(:,k);
        end
    end
end

% 3) class means (you can swap mean?median for extra robustness)
Topo.Spike    = mean(S, 2, 'omitnan');
Topo.NoSpike  = mean(N, 2, 'omitnan');
Topo.Diff     = Topo.Spike - Topo.NoSpike;    % simple contrast
Topo.refLbl   = refLbl(:);
Topo.win_ms   = win_ms;

counts = struct('Spike', size(S,2), 'NoSpike', size(N,2));
end
