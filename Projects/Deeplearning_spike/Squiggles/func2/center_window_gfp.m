function Xc = center_window_gfp(X, fs, win_ms)
W = max(1, round(win_ms*fs/1000));
g = std(X,[],1); [~,ip]=max(g); T=size(X,2);
L=max(1, ip - floor(W/2)); R=min(T, L+W-1);
Xc = X(:, L:R);
end
