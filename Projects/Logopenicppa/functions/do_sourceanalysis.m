function source_active = do_sourceanalysis(mcfg, ep_data)

cfg = [];
cfg.channel = 'MEG';
cfg.covariance = 'yes';
cov_matrix = ft_timelockanalysis(cfg, ep_data);

cfg = [];
cfg.method = 'lcmv';
cfg.grid = mcfg.individual_grid;
cfg.headmodel = mcfg.vol;
cfg.lcmv.lambda = '10%';
cfg.lcmv.keepfilter = 'yes';
source_commonfilter = ft_sourceanalysis(cfg, cov_matrix);

cfg = [];
cfg.method = 'lcmv';
cfg.grid = mcfg.individual_grid;
cfg.grid.filter = source_commonfilter.avg.filter;
cfg.headmodel = mcfg.vol;
cfg.rawtrial = 'yes';
cfg.keeptrials = 'yes';
cfg.lcmv.projectmom = 'yes';
cfg.lcmv.fixedori = 'yes';
source_active = ft_sourceanalysis(cfg, ep_data);


