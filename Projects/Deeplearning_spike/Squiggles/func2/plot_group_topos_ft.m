function plot_group_topos_ft(Topo, grad_or_sens)
% grad_or_sens: FieldTrip sensor struct (e.g., raw_data.grad or sens)
% Draws Spike mean, NoSpike mean, and Diff in one figure

% Prepare a combined layout for ALL MEG channels
cfgL = []; cfgL.grad = grad_or_sens; cfgL.channel = 'MEG';
lay  = ft_prepare_layout(cfgL);

% Make a tiny "timelock" tpl from a C×1 vector
mkTL = @(vals) struct('label',{Topo.refLbl}, 'dimord','chan_time', ...
                      'time',0,'avg',vals(:));

% Common plot cfg
base = [];
base.parameter = 'avg';
base.xlim      = [0 0];
base.layout    = lay;
base.marker    = 'off';
base.comment   = 'no';
base.colorbar  = 'yes';

figure('Color','w');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% Spike
nexttile;
cfg = base; cfg.zlim='maxabs';
ft_topoplotER(cfg, mkTL(Topo.Spike));
title(sprintf('Spike mean (N=%d, ±%d ms)', ...
    isfieldstruct(Topo,'counts','Spike'), Topo.win_ms), 'Interpreter','none');

% NoSpike
nexttile;
cfg = base; cfg.zlim='maxabs';
ft_topoplotER(cfg, mkTL(Topo.NoSpike));
title(sprintf('NoSpike mean (N=%d, ±%d ms)', ...
    isfieldstruct(Topo,'counts','NoSpike'), Topo.win_ms), 'Interpreter','none');

% Difference
nexttile;
% lock zlim to symmetric max of |diff|
cfg = base; m = max(abs(Topo.Diff)); cfg.zlim = [-m m] + 1e-12*[-1 1];
ft_topoplotER(cfg, mkTL(Topo.Diff));
title('Spike - NoSpike');
end

function out = isfieldstruct(S, f, g)
% helper to safely pull counts if you passed them; else show 0
out = 0;
if isstruct(S) && isfield(S, 'counts') && isfield(S.counts, f)
    out = S.counts.(f);
elseif nargin>=3 && isfield(S, g)
    out = S.(g);
end
end
