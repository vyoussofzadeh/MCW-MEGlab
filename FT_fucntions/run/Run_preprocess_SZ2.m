
% if exist(fullfile(outd.sub,[Datalog.spktpye,'_ic_',subj,'.mat']), 'file') == 2
%     load(fullfile(outd.sub,[Datalog.spktpye, '_ic_',subj,'.mat']));
% else   
    
%     cfg = [];
%     cfg.dataset = datafile;
%     cfg.hpfilter = 'yes';
%     cfg.hpfreq = 0.1;
%     cfg.lpfilter = 'yes';
%     cfg.lpfreq = 40;
%     cfg.channel = {'megmag', 'meggrad'};
% %     cfg.coilaccuracy = 0;
%     f_data = ft_preprocessing(cfg);
    
    cfg        = [];
    cfg.dataset = datafile;
    cfg.demean = 'yes';
    cfg.dftfilter = 'yes';
    cfg.hpfilter = 'yes';
    cfg.lpfilter = 'yes';
    cfg.hpfiltord = 3;
    cfg.lpfreq = 40;
    cfg.hpfreq = 0.1;
%     cfg.dftfreq = 50;
    cfg.channel = {'megmag', 'meggrad'};
    f_data        = ft_preprocessing(cfg);
 
    
    %%
    cfg = [];
    cfg.resamplefs = 500;
    data_resampled = ft_resampledata(cfg, f_data);

    %% ICA cleaning
    warning off
    if flag.preprocessing.ica == 1
        cfg = [];
        cfg.savepath = outd.sub;
%         cfg.savepath = fullfile(outd.sub,['ic_',subj,'.mat']);
        cfg.saveflag = 2;
        cfg.overwrite = 2;
        cfg.lay = lay;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = ic_selection;
        cfg.savefile = fullfile(cfg.savepath,[Datalog.spktpye,'_ic_',cfg.subj,'.mat']);
%         cln_data = vy_ica_cleaning_light(cfg, r_data);
        cln_data = vy_ica_cleaning_light_SZ(cfg, data_resampled);
        disp('ICA cleaning was completed');
    end
% end

    
