
if exist(fullfile(outd.sub,['Block_ic_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['Block_ic_',subj,'.mat']));
    load(fullfile(outd.sub,['Block_trial_info_updt_',subj,'.mat']));
else
    
    if flag.trialinfo.filtering == 1
        savepath = fullfile(outd.sub,['trial_info_',subj,'.mat']);
        if exist(savepath, 'file') == 2
            load(savepath);
        else
            %%
            trlInfoColDescr = [];
            
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
            save(savepath, 'trl', 'trlInfoColDescr','trialSummary','trialinfo');
        end
    end
    
    cfg                         = [];
    cfg.dataset                 = datafile;
    cfg.trialfun                = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.prestimTime        = trlInfoColDescr(:,8); % in seconds
    cfg.trialdef.poststimTime       = trlInfoColDescr(:,9); % in seconds
    cfg.trialdef.cutmode = 1;
    cfg = ft_definetrial(cfg);
    trl_sel = [unique(trlInfoColDescr(:,8)),unique(trlInfoColDescr(:,9))];    
    cfg = [];
    cfg.channel    = {'MEG'};
    cfg.dataset                 = datafile;
    data_all = ft_preprocessing(cfg);
    
    %%
    idx = find(diff(trlInfoColDescr(:,3)) == 1);
    idtask = [trlInfoColDescr(idx,2);trlInfoColDescr(end,2)];

    %%
    data_epoch1 = data_all;
    data_epoch1.trial=[];
    tinterval = 1e3;
    ttrl = data_epoch1.time{1};
    trl_res =[];
    for i = 1:length(trl_sel)
        trl_res.trial{i} = data_all.trial{1}(:, trl_sel(i,1)-tinterval:trl_sel(i,2));
        trl_res.sampleinfo(i,:) = [trl_sel(i,1)-tinterval,trl_sel(i,2)];
        trl_res.time{i} = ttrl(:, trl_sel(i,1)-tinterval:trl_sel(i,2));
    end
    trl_res.hdr = data_all.hdr;
    trl_res.label = data_all.label;
    trl_res.grad = data_all.grad;
    trl_res.fsample = data_all.fsample;
    %%
    
    %%
    if flag.preprocessing.filtering == 1
        %% Filteting, Event reading, Artifact rejecrtion
        savepath = fullfile(outd.sub,['Blk_f_',subj,'.mat']);
        if exist(savepath, 'file') == 2
            load(savepath)
        else
            %%
            disp('preprocessing ...');
            cfg = [];
            cfg.hpfilter = 'yes';
            cfg.lpfilter = 'yes';
            cfg.dftfilter = 'yes';
            cfg.hpfiltord = 3;
            cfg.hpfreq = 1;
            cfg.lpfreq = 40;
            f_data = ft_preprocessing(cfg,trl_res);
            disp('filtering was completed');
            save(savepath, 'f_data', '-v7.3');
        end
    end
        
    %% Bad trials & channels (automated)
%     if flag.preprocessing.artifact == 1
%         
%         %%
%         accuracy = length(find(trlInfoColDescr(:,14)==1))/length(trlInfoColDescr)*100;
%         disp(['Accuracy was: %', num2str(accuracy)]);
%         %%
%         cfg = [];
%         cfg.resamplefs  = 1000;
%         cfg.demean      = 'no';
%         cfg.detrend     = 'no';
%         f_data = ft_resampledata(cfg, f_data);
%         savepath = fullfile(outd.sub,['Blk_a_',subj,'.mat']);
%         cfg = [];
%         cfg.pflag = 1; % yes:1, No:2
%         cfg.saveflag = 1; % yes:1, No:2
%         cfg.savepath = savepath;
%         cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)];%[-200,900]./1000;
%         cfg.rejectpercentage = .95;
%         [r_data,report] = vy_artifactreject(cfg, f_data);
%         % disp('Bad data rejection was completed');
%     end
    
    %%
    if flag.preprocessing.ica == 1
        
        %%
%         idtask1 = idtask;
%         idtask1(report.btrl,:) = [];
        
        cfg = [];
        cfg.savepath = outd.sub;
        cfg.savefile = fullfile(outd.sub,['Blk_ic_',subj,'.mat']);
        cfg.saveflag = 1;
        cfg.overwrite = 2;
        cfg.lay = lay;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = ic_selection;
        %         cln_data = vy_ica_cleaning_light(cfg, r_data);
        [cln_data,report_ic] = vy_ica_cleaning_light3(cfg, f_data);
        disp('ICA cleaning was completed');
        
%         idtask1(report_ic.btrl,:) = [];
        savepath = fullfile(outd.sub,['Blk_trial_info_updt_',subj,'.mat']);
        save(savepath, 'idtask');
        
    end
    
    
end