function X1D = make_features_1D_spectral_flat(X, fs)
% raw + bp_beta + bp_alpha + entropy (flatten)
S = feats_spectral_epoch(X, fs);
C = size(X,1); T = size(X,2);
bp_beta_map = repmat(S.bp_beta, 1, T);
bp_alpha_map= repmat(S.bp_alpha,1, T);
sent_map    = repmat(S.sent,    1, T);
Z = [X; bp_beta_map; bp_alpha_map; sent_map];   % (4C)×T
X1D = single(Z);
end
