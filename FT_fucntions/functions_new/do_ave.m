function a_data = do_ave(data)

tt = data.time{1};
idx = tt==0;

cfg                   = [];
cfg.covariance        = 'yes';
cfg.covariancewindow  = 'all';
cfg.preproc.baselinewindow = [tt(1),tt(idx)];
cfg.preproc.demean    = 'yes';    % enable demean to remove mean value from each single trial

a_data         = ft_timelockanalysis(cfg, data);
