function ep_data = vy_epoch_paired(r_data1,r_data2, toi)

% ep_data.all = r_data;

cfg = [];
% cfg.toilim = [-0.310 -0.010];
% cfg.toilim = [-0.4 -0.1];
cfg.toilim = toi;
ep_data.bsl = ft_redefinetrial(cfg, r_data1);

% cfg.toilim = [0.5 0.8];
cfg.toilim = toi;
% cfg.toilim = [0.6 0.9];
% cfg.toilim = [0.3 0.6];
ep_data.pst = ft_redefinetrial(cfg, r_data2);

% cfg.toilim = [0 0.3];
% ep_data.aud = ft_redefinetrial(cfg, r_data);