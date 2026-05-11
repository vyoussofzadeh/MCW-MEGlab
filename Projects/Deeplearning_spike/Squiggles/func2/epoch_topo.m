function vals = epoch_topo(X, fs, win_ms, idxWin)
% X: [C x T] epoch, fs in Hz
% win_ms: scalar window width around peak GFP (e.g., 80). Ignored if idxWin provided.
% idxWin: optional logical/indices for time samples to average over (e.g., t>=0.08 & t<=0.12)

C = size(X,1);
if nargin>=4 && ~isempty(idxWin)
    vals = mean(X(:, idxWin), 2);
else
    if nargin<3 || isempty(win_ms), win_ms = 80; end
    g = std(X,[],1); [~,ip] = max(g);
    W = max(1, round(win_ms*fs/1000));
    i1 = max(1, ip - floor(W/2));
    i2 = min(size(X,2), ip + floor(W/2));
    vals = mean(X(:, i1:i2), 2);    % C×1 signed average
end
% optional: de-mean across channels to emphasize pattern
% vals = vals - median(vals,'omitnan');
end
