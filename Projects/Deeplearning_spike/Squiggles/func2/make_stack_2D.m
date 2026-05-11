function X2D = make_stack_2D(X, fs, XY, opts)
% opts: struct toggles (true/false): time, spatial, phase, longview
% returns C×T×D
if nargin<4 || isempty(opts)
    opts.time = true; opts.spatial=false; opts.phase=false; opts.longview=false;
end
C=size(X,1); T=size(X,2);
S = {}; names = {};

% raw always first
S{end+1} = X; names{end+1}='raw';

if opts.time
    M = feats_morph_time(X, fs);
    S = [S, {M.d1, M.d2, M.env}]; names = [names, {'d1','d2','env'}];
    % optional extras: M.ll, M.tke, broadcast(M.zcr)
end

if opts.spatial && ~isempty(XY)
    [Xs, Xe] = feats_spatial_graph(X, XY, 0.1);
    S = [S, {Xs, Xe}]; names = [names, {'smooth','edge'}];
end

if opts.phase && ~isempty(XY)
    P = feats_phase_sync(X, fs, 3, XY);
    S{end+1} = P.lcv; names{end+1}='lcv';
end

% concatenate depth
X2D = single(cat(3, S{:}));
end

function X1D = make_features_1D_from_stack(X2D)
% C×T×D -> (C*D)×T
X1D = reshape(X2D, size(X2D,1)*size(X2D,3), size(X2D,2));
end
