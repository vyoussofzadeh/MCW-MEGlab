function ep_data = do_epoch(r_data,toi)

ep_data.all = r_data;

cfg = [];
cfg.toilim = toi(1,:);
ep_data.bsl = ft_redefinetrial(cfg, r_data);

cfg.toilim = toi(2,:);
ep_data.pst = ft_redefinetrial(cfg, r_data);