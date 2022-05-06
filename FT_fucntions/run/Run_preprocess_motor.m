
if exist(fullfile(outd.sub,['icl_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['icl_',subj,'.mat'])); cln.left = cln_data;
    load(fullfile(outd.sub,['icr_',subj,'.mat'])); cln.right = cln_data;    
else
    
    if flag.preprocessing.filtering == 1
        %% Filteting, Event reading, Artifact rejecrtion
        savepath1 = fullfile(outd.sub,['fl_',subj,'.mat']);
        savepath2 = fullfile(outd.sub,['fr_',subj,'.mat']);
        if exist(savepath1, 'file') == 2
            load(savepath1)
            load(savepath2)
        else
%             cfg = [];
%             cfg.datafile  = datafile;
%             cfg.hpfreq = 0.1;
%             cfg.lpfreq = 40;
%             cfg.hpfiltord = 3;
%             cfg.hpfilter = 'yes';
%             cfg.lpfilter = 'yes';
%             cfg.dftfilter = 'yes';
%             [f_data] = ft_preprocessing(cfg);
%             
%             cfg                         = [];
%             cfg.datafile  = datafile;
%             cfg.channel                 = {'MEG'};
%             cfg.baselinewindow          = [-0.45 0.0];
%             cfg.trialfun                = 'ft_trialfun_general'; % this is the default
%             cfg.trialdef.eventtype      = epoch_type;
%             cfg.trialdef.eventvalue     = 16; % the value of the stimulus trigger for fully incongruent (FIC).
%             cfg.trialdef.prestim        = 1; % in seconds
%             cfg.trialdef.poststim       = 3; % in seconds
%             cfg = ft_definetrial(cfg);
%             trl_left = cfg.trl;
%             cfg.trialdef.eventvalue     = 32; % the value of the stimulus trigger for fully incongruent (FIC).
%             cfg = ft_definetrial(cfg);
%             trl_right = cfg.trl;
%             
%             cfg = [];
%             cfg.trl = trl_left;
%             fl_data = ft_redefinetrial(cfg,f_data);
%             cfg.trl = trl_right;
%             fr_data = ft_redefinetrial(cfg,f_data);
            
            
            cfg = [];
            cfg.eventid = 16;
            cfg.epochtype = epoch_type;
            cfg.datafile  = datafile;
            cfg.hpfreq = 0.1;
            cfg.lpfreq = 40;
            [fl_data] = vy_preprocess(cfg);
            cfg.eventid = 32;
            [fr_data] = vy_preprocess(cfg);
            disp('filtering was completed');
            save(savepath1, 'fl_data', '-v7.3');
            save(savepath2, 'fr_data', '-v7.3');
        end
    end
    
    %%
    %     cfg = [];
    %     cfg.resamplefs  = 250;
    %     cfg.demean      = 'no';
    %     cfg.detrend     = 'no';
    %     f_data = ft_resampledata(cfg, f_data);
    
    %% Bad trials & channels (automated)
    if flag.preprocessing.artifact == 1
        
        savepath1 = fullfile(outd.sub,['al_',subj,'.mat']);
        savepath2 = fullfile(outd.sub,['ar_',subj,'.mat']);
        cfg = [];
        cfg.pflag = 2; % yes:1, No:2
        cfg.saveflag = 2; % yes:1, No:2
        cfg.savepath = savepath1;
        cfg.latency = [fl_data.time{1}(1),fl_data.time{1}(end)];%[-200,900]./1000;
        cfg.rejectpercentage = .95;
        [rl_data,report_l] = vy_artifactreject(cfg, fl_data);
        cfg.latency = [fr_data.time{1}(1),fr_data.time{1}(end)];%[-200,900]./1000;
        cfg.savepath = savepath2;
        [rr_data,report_r] = vy_artifactreject(cfg, fr_data);
        % disp('Bad data rejection was completed');
    end
    %% Bad trials & channels (Manuual)
    % clear r_data
    % cfg = [];
    % cfg.metric = 'zvalue';  % use by default zvalue method
    % cfg.latency = [-200,900];
    % cfg.layout   = lay;   % this allows for plotting individual trials
    % r_data   = ft_rejectvisual(cfg, f_data);
    
    %% Inspecting bad data
    % cfg = [];
    % cfg.viewmode = 'vertical';
    % cfg.continuous = 'no';
    % cfg.trials     = report.btrl;
    % cfg.channel   = report.bchan;
    % ft_databrowser(cfg,f_data);
    
    %% ICA cleaning
    if flag.preprocessing.ica == 1
        cfg = [];
        cfg.savepath = outd.sub;
        %         cfg.savepath = fullfile(outd.sub,['ic_',subj,'.mat']);
        cfg.saveflag = 1;
        cfg.overwrite = 2;
        cfg.lay = lay;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = ic_selection;
        cfg.savefile = fullfile(cfg.savepath,['icl_',cfg.subj,'.mat']);
        %         cln_data = vy_ica_cleaning_light(cfg, r_data);
        clnl_data = vy_ica_cleaning_light2(cfg, rl_data);
        cfg.savefile = fullfile(cfg.savepath,['icr_',cfg.subj,'.mat']);
        clnr_data = vy_ica_cleaning_light2(cfg, rr_data);
        disp('ICA cleaning was completed');
    end
end