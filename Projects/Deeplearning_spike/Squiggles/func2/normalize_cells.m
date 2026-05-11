function Xout = normalize_cells(Xin, stats, mode)
% NORMALIZE_CELLS  Per-channel normalization for {C×T} epochs.
%   Xout = normalize_cells(Xin, stats, 'zscore')
%   Xout = normalize_cells(Xin, stats, 'robust')
%
% Inputs
%   Xin  : 1×N cell, each C×T numeric
%   stats: struct from fit on TRAIN only:
%          - zscore: stats.mu (C×1), stats.sd (C×1)
%          - robust: stats.med (C×1), stats.madv (C×1)
%   mode : 'zscore' or 'robust'
%
% Output
%   Xout : 1×N cell, same sizes as Xin

assert(iscell(Xin) && ~isempty(Xin), 'Xin must be a non-empty cell array.');
C = size(Xin{1},1);
switch lower(mode)
    case 'zscore'
        mu = stats.mu(:);  sd = stats.sd(:);
        assert(numel(mu)==C && numel(sd)==C, 'stats size mismatch.');
        f = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mu), max(sd,1e-12));
    case 'robust'
        med = stats.med(:);  madv = stats.madv(:);
        assert(numel(med)==C && numel(madv)==C, 'stats size mismatch.');
        f = @(x) bsxfun(@rdivide, bsxfun(@minus, x, med), max(madv,1e-12));
    otherwise
        error('normalize_cells:UnknownMode','Use ''zscore'' or ''robust''.');
end

Xout = cell(size(Xin));
for i = 1:numel(Xin)
    Xi = Xin{i};
    if isempty(Xi), Xout{i} = Xi; continue; end
    if size(Xi,1) ~= C, error('Epoch %d has C=%d (expected %d).', i, size(Xi,1), C); end
    Xout{i} = f(Xi);
end
end
