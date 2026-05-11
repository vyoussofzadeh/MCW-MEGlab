function [Z, mu, sig] = sanitize_features(X, mu, sig)
% SANITIZE_FEATURES  Replace NaN/Inf and z-score columns.
% Usage:
%   [Fz, mu, sig] = sanitize_features(Ftr);   % fit on TRAIN
%   Fvaz         = sanitize_features(Fva, mu, sig);  % apply to VAL/TEST

% 1) clean
X(~isfinite(X)) = 0;

% 2) fit stats on TRAIN if not provided
if nargin < 2 || isempty(mu)
    mu = mean(X, 1, 'omitnan');
    sig = std(X, 0, 1, 'omitnan');
    sig(sig < 1e-8) = 1;
end

% 3) standardize
Z = (X - mu) ./ sig;
end
