
if exist(fullfile(outd.sub,['ic_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['ic_',subj,'.mat']));
    load(fullfile(outd.sub,['cond_',subj,'.mat']));
else
    
    
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
    cfg.demean = 'yes';
    cfg.baselinewindow = [-0.45 0.0];
    f_data = ft_preprocessing(cfg);

    %% ICA
    if flag.preprocessing.ica == 1
        cfg = [];
        cfg.savepath = outd.sub;
        cfg.savefile = fullfile(outd.sub,['ic_',subj,'.mat']);
        cfg.saveflag = 1;
        cfg.overwrite = 1;
        cfg.lay = lay;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = ic_selection;
        cln_data = vy_ica_cleaning_light(cfg, f_data);
        disp('ICA cleaning was completed');
    end
    
    %%
    savepath = fullfile(outd.sub,['cond_',subj,'.mat']);
    save(savepath, 'trl');
    
end

