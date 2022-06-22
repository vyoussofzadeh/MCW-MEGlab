
%- Cov matrix
cfg = [];
cfg.channel = 'MEG';
cfg.covariance = 'yes';
cov_matrix = ft_timelockanalysis(cfg, ep_data.all);

cfg = [];
cfg.method = 'lcmv';
cfg.grid = individual_grid;
cfg.headmodel = individual_headmodel;
cfg.lcmv.lambda = '10%';
cfg.lcmv.keepfilter = 'yes';
source_whole = ft_sourceanalysis(cfg, cov_matrix);

cfg = [];
cfg.method = 'lcmv';
cfg.grid = individual_grid;
cfg.grid.filter = source_whole.avg.filter;
cfg.headmodel = individual_headmodel;
cfg.rawtrial = 'yes';
cfg.keeptrials = 'yes';
cfg.lcmv.projectmom = 'yes';
cfg.lcmv.fixedori = 'yes';
source_active = ft_sourceanalysis(cfg, ep_data.pst);
