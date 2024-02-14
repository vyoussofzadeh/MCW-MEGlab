function f_data = do_ft_preprocess(cfg_in)

cfg                         = [];
cfg.dataset                 = cfg_in.datafile;
cfg.trialfun                = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype      = 'STI101';
cfg.trialdef.eventvalue     = cfg_in.eventvalue; % the value of the stimulus trigger for fully incongruent (FIC).
cfg.trialdef.prestim        = cfg_in.prestim; % in seconds
cfg.trialdef.poststim       = cfg_in.poststim; % in seconds
cfg = ft_definetrial(cfg);

cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.hpfiltord = 3;
cfg.hpfreq = cfg_in.hpfreq;
cfg.lpfreq = cfg_in.lpfreq;
cfg.channel = {'MEG'};
cfg.dftfreq = cfg_in.dftfreq;
cfg.dftfilter = 'yes';
f_data = ft_preprocessing(cfg);


end