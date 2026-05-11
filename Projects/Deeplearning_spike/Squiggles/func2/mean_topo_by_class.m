function M = mean_topo_by_class(Xcells, Y, fs, refLbl, win_ms)
% Returns struct with mean topographies per class and per type (mag, grad)
cats = categories(Y);
isMag  = endsWith(refLbl,'1');
isGrad = endsWith(refLbl,'2') | endsWith(refLbl,'3');

% Collect per-epoch topographies
Tmag = []; Tgrad = []; L = [];  % rows = channels, cols = epochs
for i = 1:numel(Xcells)
    v = epoch_topo(Xcells{i}, fs, win_ms, []);
    Tmag  = [Tmag,  v(isMag )]; %#ok<AGROW>
    Tgrad = [Tgrad, v(isGrad)]; %#ok<AGROW>
    L     = [L;     Y(i)];      %#ok<AGROW>
end

% Means per class
M = struct();
for c = 1:numel(cats)
    msk = (L == cats{c});
    M.(cats{c}).mag  = mean(Tmag(:, msk),  2, 'omitnan');   % [Nmag x 1]
    M.(cats{c}).grad = mean(Tgrad(:, msk), 2, 'omitnan');   % [Ngrad x 1] (raw separate coils)
end
M.isMag = isMag; M.isGrad = isGrad; M.refLbl = refLbl; M.win_ms = win_ms; M.cats = cats;
end
