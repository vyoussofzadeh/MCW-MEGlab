function X2D = make_stack_2D_spectral_flat(X, fs)
% X: C×T -> C×T×D where D= (raw + beta + alpha + entropy)
S = feats_spectral_epoch(X, fs);
M = spectral_maps_from_scalars(S, size(X,2)); % all C×T
% choose a compact set (you can add delta/theta/gamma if useful)
X2D = single(cat(3, X, M.bp_beta, M.bp_alpha, M.sent));
% depth order: {'raw','bp_beta','bp_alpha','sent'}
end
