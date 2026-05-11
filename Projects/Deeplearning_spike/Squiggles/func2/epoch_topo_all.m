function vals = epoch_topo_all(X, fs, win_ms)
% X: [C x T] one epoch; fs: Hz; win_ms: width around GFP peak (e.g., 80)
if nargin<3 || isempty(win_ms), win_ms = 80; end
g  = std(X,[],1); [~,ip] = max(g);
W  = max(1, round(win_ms*fs/1000));
i1 = max(1, ip - floor(W/2));
i2 = min(size(X,2), ip + floor(W/2));
vals = mean(X(:, i1:i2), 2);          % Cx1 signed mean (use abs(...) for magnitude map)
% Optional: center across channels to emphasize pattern
% vals = vals - median(vals,'omitnan');
end
