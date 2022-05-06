
if exist(fullfile(outd.sub,['ic_str_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['ic_str_',subj,'.mat']));
    load(fullfile(outd.sub,['ic_math_',subj,'.mat']));
else
    
    if flag.preprocessing.filtering == 1
        %% Filteting, Event reading, Artifact rejecrtion
        savepath = fullfile(outd.sub,['f_',subj,'.mat']);
        if exist(savepath, 'file') == 2
            load(savepath)
%             load('trial_info')
        else
            %             event = ft_read_event(datafile);
            %             clear val
            %             for i=1:length(event)
            %                 val(i) = event(i).value;
            %             end
            %             val1 = unique(val);
            %             Evnt_IDs = mode(val);
            %
            %
            %             a = unique(val);
            %             n = histc(val,a);
            %             [n,idx] = sort(n);
            %             myFreq = a(idx); %Unique Values
            
            %%
            scanmnem = 'StoryM_neuromag';
            trialFunName=['trialfun_',scanmnem];
            trialDefFuncFile=['alltrialdefparams_',scanmnem];
            trialDefFuncHandle=str2func(trialDefFuncFile);
            savesuffix_trialinfo='trialinfo';
            
            experimentid = 'SM';
            scanid = run;
            pipelineprefix='baddata';
            
            resultprefix = sprintf('%s_%s_%s', experimentid, scanid, pipelineprefix);
            %------------------------------------------
            %--- Extract trials definition
            allTrlCfgs=trialDefFuncHandle();
            Ncasegroups=length(allTrlCfgs);
            outputTrialInfoSummaryFile = [resultprefix,'_raw',savesuffix_trialinfo,'_QC.txt'];
            montage         = hcp_exgmontage(subj, experimentid, scanid);
            iCaseGroup=1;
            
            %------------------------------------------
            trlCfg                 = allTrlCfgs{iCaseGroup};
            trlCfg.datafile        = datafile;
            trlCfg.trialfun        = trialFunName;
            if ~isempty(regexp(scanid,'Motor', 'once'))
                trlCfg.trialdef.montage=montage;
            end
            %------------------------------------------
            eval(['[trl,trialinfo,trlInfoColDescr,trialSummary,scanStartSamp,scanEndSamp,warninfo]=',trialFunName,'(trlCfg);']);
            % Save the summary ascii file
            hcp_write_ascii(outputTrialInfoSummaryFile,'trialSummary'); % Trial Summary should be the same no matter what the input arguments are. This is because the trial definition function creates all information about the trial definition. This is what the summary contains. The trl field only contains the trials defined by the input arguments.
            ErrorFile=['ERROR_',outputTrialInfoSummaryFile];
            delete([ErrorFile,'*']); % Remove any error files present from previous runs
            WarningFile=['WARNING_',outputTrialInfoSummaryFile];
            if isempty(warninfo)
                delete([WarningFile,'*']); % Remove any warning files present from previous runs
            else
                hcp_write_ascii(WarningFile,'warninfo');
            end
            save('trial_info', 'trl', 'trlInfoColDescr','trialSummary','trialinfo');
            
            %%
            disp('preprocessing ...');
            cfg                         = [];
            cfg.dataset                 = datafile;
            cfg.trialfun                = 'trialfun_StoryM_BaseExtractAll_neuromag'; % this is the default
            cfg.trialdef.prestimTime        = 1; % in seconds
            cfg.trialdef.poststimTime       = 3; % in seconds
            cfg.trialdef.cutmode = 1;
            cfg = ft_definetrial(cfg);
            
            cfg.hpfilter = 'yes';
            cfg.lpfilter = 'yes';
            cfg.dftfilter = 'yes';
            cfg.hpfiltord = 3;
            % cfg.dftfreq = [60 120 180];
            cfg.hpfreq = 1;
            cfg.lpfreq = 40;
            % cfg.hpfreq = 2;
            % cfg.lpfreq = 30;
            cfg.channel = {'MEG'};
            % cfg.channel = {'MEGGRAD'};
            cfg.demean = 'yes';
            cfg.baselinewindow = [-0.45 0.0];
            f_data = ft_preprocessing(cfg);
            disp('filtering was completed');
            save(savepath, 'f_data', '-v7.3');
        end
    end
    
    %%
    accuracy = length(find(trlInfoColDescr(:,14)==1))/length(trlInfoColDescr)*100;
    disp(['Accuracy was: %', num2str(accuracy)])
    
    %%
    idx_str = find(trlInfoColDescr(:,2)==1);
    idx_math = find(trlInfoColDescr(:,2)==2);
    
    cfg = [];
    cfg.trials = idx_str;
    f_data_str = ft_selectdata(cfg,f_data);
    cfg.trials = idx_math;
    f_data_math = ft_selectdata(cfg,f_data);
    
    %%
    cfg = [];
    cfg.resamplefs  = 1000;
    cfg.demean      = 'no';
    cfg.detrend     = 'no';
    f_data_math = ft_resampledata(cfg, f_data_math);
    f_data_str = ft_resampledata(cfg, f_data_str);
    
    %%
    cfg = [];
    cfg.toilim = [-.4,2];
    f_data_str = ft_redefinetrial(cfg, f_data_str);
    
    cfg = [];
    cfg.toilim = [-.4,2];
    f_data_math = ft_redefinetrial(cfg, f_data_math);
    
    %% Bad trials & channels (automated)
    if flag.preprocessing.artifact == 1
        savepath = fullfile(outd.sub,['a_str_',subj,'.mat']);
        data_in = f_data_str;
        cfg = [];
        cfg.pflag = 2; % yes:1, No:2
        cfg.saveflag = 1; % yes:1, No:2
        cfg.savepath = savepath;
        cfg.latency = [data_in.time{1}(1),data_in.time{1}(end)];%[-200,900]./1000;
        cfg.rejectpercentage = .95;
        [r_data_str,report] = vy_artifactreject(cfg, data_in);
        % disp('Bad data rejection was completed');
    end
    
    if flag.preprocessing.artifact == 1
        savepath = fullfile(outd.sub,['a_math_',subj,'.mat']);
        data_in = f_data_math;
        cfg = [];
        cfg.pflag = 2; % yes:1, No:2
        cfg.saveflag = 1; % yes:1, No:2
        cfg.savepath = savepath;
        cfg.latency = [data_in.time{1}(1),data_in.time{1}(end)];%[-200,900]./1000;
        cfg.rejectpercentage = .95;
        [r_data_math,report] = vy_artifactreject(cfg, data_in);
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
    %     cfg = [];
    %     cfg.viewmode = 'vertical';
    %     cfg.continuous = 'no';
    %     cfg.trials     = report.btrl;
    %     cfg.channel   = report.bchan;
    %     ft_databrowser(cfg,f_data);
    %
    %% ICA cleaning
    if flag.preprocessing.ica == 1
        cfg = [];
        cfg.savepath = outd.sub;
        cfg.savefile = fullfile(outd.sub,['ic_str_',subj,'.mat']);
        cfg.saveflag = 1;
        cfg.overwrite = 2;
        cfg.lay = lay;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = ic_selection;
        %         cln_data = vy_ica_cleaning_light(cfg, r_data);
        cln_data_str = vy_ica_cleaning_light2(cfg, r_data_str);
        disp('ICA cleaning was completed');
    end
    
    if flag.preprocessing.ica == 1
        cfg = [];
        cfg.savepath = outd.sub;
        cfg.savefile = fullfile(outd.sub,['ic_math_',subj,'.mat']);
        cfg.saveflag = 1;
        cfg.overwrite = 2;
        cfg.lay = lay;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = ic_selection;
        %         cln_data = vy_ica_cleaning_light(cfg, r_data);
        cln_data_math = vy_ica_cleaning_light2(cfg, r_data_math);
        disp('ICA cleaning was completed');
    end
    
    
    
end