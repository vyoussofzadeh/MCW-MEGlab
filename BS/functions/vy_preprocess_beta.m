function f_data = vy_preprocess_beta(datafile,Evnt_IDs,epoch_type)


cfg                         = [];
cfg.dataset                 = datafile;
cfg.trialfun                = 'ft_trialfun_general'; % this is the default
cfg.trialdef.eventtype      = epoch_type;
cfg.trialdef.eventvalue     = Evnt_IDs; % the value of the stimulus trigger for fully incongruent (FIC).
cfg.trialdef.prestim        = 1; % in seconds
cfg.trialdef.poststim       = 2; % in seconds
cfg = ft_definetrial(cfg);

cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.dftfilter = 'yes';
cfg.hpfiltord = 3;
% cfg.dftfreq = [60 120 180];
cfg.hpfreq = 12;
cfg.lpfreq = 25;
% cfg.hpfreq = 2;
% cfg.lpfreq = 30;
% cfg.bsfreq = 32;
cfg.channel = {'MEG'};
% cfg.channel = {'MEGGRAD'};
cfg.demean = 'yes';
cfg.baselinewindow = [-0.45 0.0];
f_data = ft_preprocessing(cfg);

% [cfg, artifact] = ft_artifact_jump(cfg);        % id jumps
% f_data = ft_rejectartifact(cfg, f_data);        % remove scanner (MEG) jumps

% cfg = [];
% cfg.dataset  = datafile;
% cfg.hpfilter = 'yes';
% cfg.lpfilter = 'yes';
% cfg.dftfilter = 'yes';
% cfg.hpfiltord = 3;
% cfg.dftfreq = [60 120 180];
% cfg.hpfreq = .1;
% cfg.lpfreq = 40;
% cfg.channel = {'MEG'};
% cfg.demean = 'yes';
% % cfg.baselinewindow = [-0.2 0.0];
% data = ft_preprocessing(cfg);
% 
% cfg = [];
% cfg.dataset                 = datafile;
% cfg.trialfun                = 'ft_trialfun_general'; % this is the default
% cfg.trialdef.eventtype      = epoch_type;
% cfg.trialdef.eventvalue     = Evnt_IDs; % the value of the stimulus trigger for fully incongruent (FIC).
% cfg.trialdef.prestim        = 1; % in seconds
% cfg.trialdef.poststim       = 2; % in seconds
% cfg = ft_definetrial(cfg);
% 
% trl = cfg.trl;
% cfg = [];
% cfg.trl = trl;
% f_data = ft_redefinetrial(cfg,data);
% 
% 
% cfg = [];
% cfg.dataset = datafile;
% cfg.trl = trl;
% [cfg, artifact] = ft_artifact_jump(cfg);        % id jumps
% f_data = ft_rejectartifact(cfg, f_data);        % remove scanner (MEG) jumps

