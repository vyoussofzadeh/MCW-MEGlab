
if exist(savepath, 'file') == 2
    load(savepath);
else
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
    
    %% Artifact rejecrtion:  removing bad trials & sensors (semi-automated)
    savepath = fullfile(outd.sub,['a_',subj,'.mat']);
    cfg = [];
    cfg.pflag = 2; % Yes:1, No:2
    cfg.saveflag = 2; % Yes:1, No:2
    cfg.savepath = savepath;
    cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)];
    cfg.rejectpercentage = .95;
    cfg.rbadtrl = 2;
    cfg.rbadsen = 1;
    r_data = vy_artifactreject2(cfg, f_data);
    
    %% ICA
    if flag.preprocessing.ica == 1
        cfg = [];
        cfg.savepath = outd.sub;
        cfg.savefile = savepath;
        cfg.saveflag = 1;
        cfg.overwrite = 1;
        cfg.lay = lay;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = ic_selection;
        cln_data = vy_ica_cleaning_light(cfg, r_data);
        disp('ICA cleaning was completed');
    end
    %     save(savepath, 'cln_data');
end

