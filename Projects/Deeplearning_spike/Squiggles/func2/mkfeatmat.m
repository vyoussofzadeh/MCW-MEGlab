function [F, keep] = mkfeatmat(Xcells, T)
% Build N×6 feature matrix from a cell of epochs; skip failures.
n = numel(Xcells);
F = zeros(0,6); keep = false(n,1);
for i = 1:n
    try
        f = feat_vs_template(single(Xcells{i}), T);   % 1×6
        if isempty(f), continue; end
        f = reshape(double(f),1,[]);                 % force 1×6
        f(~isfinite(f)) = 0;                         % clean NaN/Inf
        F(end+1, :) = f;                             %#ok<AGROW>
        keep(i) = true;
    catch
        % skip this epoch
    end
end
end
