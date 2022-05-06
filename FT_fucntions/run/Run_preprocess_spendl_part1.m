
%% Filteting, Event reading
disp('preprocessing ...');
cfg                         = [];
cfg.dataset                 = datafile;
cfg.trialfun                = 'trialfun_sentence_meaningfulness_SP'; % this is the default
cfg.oneset_type = oneset_type; % Choices: 'sentence', 'fixation', 'target', 'response'
cfg.plotflag = 2;
cfg.trialdef.prestimTime        = prestimTime; % in seconds
cfg.trialdef.poststimTime       = poststimTime; % in seconds
%         cfg.trials = 1:5;
cfg = ft_definetrial(cfg);

idx  = find(cfg.trl(:,1) < 0);
if ~isempty(idx)
    cfg.trl(idx,:) = [];
end

trl = cfg;

cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.dftfilter = 'yes';
cfg.hpfiltord = 3;
cfg.hpfreq = 1;
cfg.lpfreq = 40;
cfg.channel = {'MEG'};
% cfg.demean = 'yes';
% cfg.baselinewindow = [-0.45 0.0];
f_data = ft_preprocessing(cfg);

%%
if f_data.fsample == 2e3
    cfg_resamp = [];
    cfg_resamp.resamplefs = 1000;
    f_data = ft_resampledata(cfg_resamp, f_data);
end
