function ft_plot_mean_topo_mag(meanVals, refLbl, titleStr)
tl = []; tl.label = refLbl(:); tl.dimord='chan_time'; tl.time=0;
tl.avg = zeros(1, numel(refLbl));  % init
isMag  = endsWith(refLbl,'1');
tl.avg(1,isMag) = meanVals(:)';    % place only mag channels

cfg = [];
cfg.parameter = 'avg'; cfg.xlim=[0 0];
cfg.layout = 'neuromag306mag.lay';
cfg.marker='off'; cfg.comment='no'; cfg.colorbar='yes'; cfg.zlim='maxabs';
ft_topoplotER(cfg, tl); title(titleStr);
end
