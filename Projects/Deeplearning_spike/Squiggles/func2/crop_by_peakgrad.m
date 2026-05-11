function Xcut = crop_by_peakgrad(X, win)
% Crop each {C×T} epoch around the strongest temporal change.
% win: target samples to keep (e.g., 120)
Xcut = cell(size(X));
for i=1:numel(X)
    xi = X{i}; if isempty(xi), Xcut{i}=xi; continue; end
    d = sum(abs(diff(xi,1,2)),1);              % edge energy across channels
    [~,t0] = max(d);
    s = max(1, t0 - floor(win/2));
    e = min(size(xi,2), s + win - 1);
    Xcut{i} = xi(:, s:e);
end
end
