function tfr = do_tfr(cfg_main, data)

% do tfr-decomposition
cfg = [];
cfg.output     = 'pow';
cfg.channel    = 'all';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
% cfg.taper      = 'dpss';
cfg.foi        = 1:3:40;
cfg.keeptrials = 'yes';
cfg.t_ftimwin  = 3./cfg.foi;
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.tapsmofrq  = 0.8 *cfg.foi;
% cfg.toi        = -0.5:0.05:3;
cfg.toi        = cfg_main.toi(1):0.05:cfg_main.toi(2);
cfg.pad = 'nextpow2';
tfr        = ft_freqanalysis(cfg, data);
set(gcf,'name',cfg_main.subj,'numbertitle','off')
if isempty(cfg_main.savefile) == 0
    save(cfg_main.savefile, 'tfr', '-v7.3');
end