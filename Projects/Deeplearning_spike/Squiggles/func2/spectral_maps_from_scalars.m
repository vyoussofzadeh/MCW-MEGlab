function M = spectral_maps_from_scalars(S, T)
% broadcast scalars to maps (C×T) so you can stack as depth for 2D CNN
fn = fieldnames(S);
for i=1:numel(fn)
    v = S.(fn{i});                        % C×1
    if isvector(v), M.(fn{i}) = repmat(v, 1, T);  % C×T
    else,            M.(fn{i}) = v;
    end
end
end
