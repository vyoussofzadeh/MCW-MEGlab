if exist(fullfile(outd.sub,['cond_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['cond_',subj,'.mat']));
    load(fullfile(outd.sub,['ic_',subj,'.mat']));
else
    
    if flag.preprocessing.filtering == 1
        %% Filteting, Event reading, Artifact rejecrtion
        savepath = fullfile(outd.sub,['f_',subj,'.mat']);
        if exist(savepath, 'file') == 2
            load(savepath)
        else
            
            event = ft_read_event(datafile,'type',{'STI101'},'value',[1,2]);
            clear val
            for i=1:length(event)
                val(i) = event(i).value;
            end
            val1 = unique(val);
            Evnt_IDs = mode(val);

            disp('preprocessing ...');
            cfg = [];
            cfg.eventid = val1;
            cfg.epochtype = event(20).type;
            cfg.datafile  = datafile;
            cfg.hpfreq = 0.1;
            cfg.lpfreq = 70;
            [f_data] = vy_preprocess(cfg);
            disp('filtering was completed');
            save(savepath, 'f_data', '-v7.3');
        end
    end    
    %%
    cfg = [];
    cfg.resamplefs  = 1000;
    cfg.demean      = 'no';
    cfg.detrend     = 'no';
    f_data = ft_resampledata(cfg, f_data);
      
    %% Bad trials & channels (automated)
    if flag.preprocessing.artifact == 1
        savepath = fullfile(outd.sub,['a_',subj,'.mat']);
        cfg = [];
        cfg.pflag = 2; % yes:1, No:2
        cfg.saveflag = 1; % yes:1, No:2
        cfg.savepath = savepath;
        cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)];%[-200,900]./1000;
        cfg.rejectpercentage = .95;
        [r_data,report] = vy_artifactreject(cfg, f_data);
        % disp('Bad data rejection was completed');
    end
    %%
    val2 = val;
    val2(report.btrl) = [];
%     
    %% ICA cleaning
    if flag.preprocessing.ica == 1
        cfg = [];
        cfg.savepath = outd.sub;
        cfg.savefile = fullfile(outd.sub,['ic_',subj,'.mat']);
        cfg.saveflag = 1;
        cfg.overwrite = 2;
        cfg.lay = lay;
        cfg.run = run;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = ic_selection;
        %         cln_data = vy_ica_cleaning_light(cfg, r_data);
        [cln_data,report_ic] = vy_ica_cleaning_light2(cfg, r_data);
%         cln_data = vy_ica_cleaning_light6(cfg, r_data);
        disp('ICA cleaning was completed');
    end
    val2(report_ic.btrl) = [];
    savepath = fullfile(outd.sub,['cond_',subj,'.mat']);
    save(savepath, 'event','val1','val2');
end