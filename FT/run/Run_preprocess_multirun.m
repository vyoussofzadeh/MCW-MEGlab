
if exist(fullfile(outd.sub,['ic_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['ic_',subj,'.mat']));
else
    
    for k = 1:nrun
        
        datafile = datafile2{k};
        
        %% Filtering
        if flag.preprocessing.filtering == 1
            %% Filteting, Event reading, Artifact rejecrtion
            savepath = fullfile(outd.sub,[num2str(k),'_f_',subj,'.mat']);
            if exist(savepath, 'file') == 2
                load(savepath)
            else
                event = ft_read_event(datafile);
                clear val
                for i=1:length(event)
                    val(i) = event(i).value;
                end
                val1 = unique(val);
                Evnt_IDs = mode(val);
                disp('preprocessing ...');
                cfg = [];
                cfg.eventid = Evnt_IDs;
                cfg.epochtype = event(20).type;
                cfg.datafile  = datafile;
                cfg.hpfreq = 0.1;
                cfg.lpfreq = 40;
                [f_data, ecg_data] = vy_preprocess(cfg);
                disp('filtering was completed');
                save(savepath, 'f_data', '-v7.3');
            end
        end
        f_data_run{k} = f_data;
    end
    f_data =  ft_appenddata([], f_data_run{:});
    f_data.grad = f_data_run{1}.grad;
    
    %% Bad trials & channels (automated)
    if flag.preprocessing.artifact == 1
        savepath = fullfile(outd.sub,['a_',subj,'.mat']);
        cfg = [];
        cfg.pflag = 2; % yes:1, No:2
        cfg.saveflag = 1; % yes:1, No:2
        cfg.savepath = savepath;
        cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)];%[-200,900]./1000;
        [r_data,report] = vy_artifactreject(cfg, f_data);
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
        cfg.savepath = fullfile(outd.sub,['ic_',subj,'.mat']);
        cfg.saveflag = 1;
        cfg.lay = lay;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        %         cln_data = vy_ica_cleaning_light(cfg, r_data);
        cln_data = vy_ica_cleaning_light2(cfg, r_data);
        disp('ICA cleaning was completed');
    end
end
