
%% Filtering, Event reading
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
% cfg.bsfilter = 'yes';
cfg.hpfiltord = 3;
cfg.hpfreq = 1;
cfg.lpfreq = 250;
% fsb = input('Enter the sop-band frequency?');
% fsb = 60
% cfg.bsfreq = [fsb-1 fsb+1]; % or whatever you deem appropriate
cfg.dftfreq = [60 120 180];
cfg.channel = {'MEG'};
% cfg.demean = 'yes';
% cfg.baselinewindow = [-0.45 0.0];
f_data = ft_preprocessing(cfg);

%sanity  check
% cfg = [];
% cfg.savefile = [];
% cfg.saveflag = 2;
% cfg.foilim = [1 250];
% cfg.plotflag  = 1;
% cfg.tapsmofrq = 5;
% cfg.taper     = 'hanning';
% vy_fft(cfg, f_data);


%%
if f_data.fsample == 2e3
    cfg_resamp = [];
    cfg_resamp.resamplefs = 1000;
    f_data = ft_resampledata(cfg_resamp, f_data);
end


