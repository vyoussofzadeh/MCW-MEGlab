
if exist(fullfile(outd.sub,['ic_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['ic_',subj,'.mat']));
%     load(fullfile(outd.sub,['cond_',subj,'.mat']));
else
    event = ft_read_event(datafile);
    clear val
    for i=1:length(event)
        val(i) = event(i).value;
        samp(i)= event(i).sample;
    end
    %             val = val - min(val);
    val1 = unique(val);
    Evnt_IDs = mode(val);
    [M,F,C] = mode(val);
%     figure,plot(val)
    
    
    if flag.preprocessing.filtering == 1
        %% Filteting, Event reading, Artifact rejecrtion
        savepath = fullfile(outd.sub,['f_',subj,'.mat']);
        if exist(savepath, 'file') == 2
            load(savepath)
        else
            
            %%
            hdr = ft_read_header(datafile);
            Fsample = hdr.Fs;
            
            Index = strfind(hdr.label,{'STI101'});
            Index = find(not(cellfun('isempty',Index)));
            
            detTrig=ft_read_data(datafile,'chanindx',Index,'header',hdr,'eventformat','neuromag_fif','dataformat','neuromag_fif');
            detTrig = (detTrig - min(detTrig));
            detTrig=bitand(detTrig,255);
            figure,plot(detTrig)
            
            detResp=ft_read_data(datafile,'chanindx',Index,'header',hdr,'eventformat','digital trigger','dataformat','neuromag_fif');
            detResp = detResp - detTrig;
            detResp = (detResp - min(detResp));
            figure,plot(detResp)
            Nsamps=length(detTrig);
            
            figure,plot(detTrig)
            for i=1:length(sampletrig)
                hold on
                y = ylim; % current y-axis limits
                plot([sampletrig(i) sampletrig(i)],y)
            end
            saveas(gcf, 'trigger_time.png')
            disp('check triggers, then enter to proceed!')
            pause,
            
            %% response/reaction time
%             resp = (detResp - mean(detResp))./max(detTrig); 
%             difResp=diff(resp);
%             difResp = abs(difResp);
%             
%             figure,plot(detTrig)
%             hold on
%             plot(resp, 'r')
%             for i=1:length(sampletrig)
%                 hold on
%                 y = ylim; % current y-axis limits
%                 plot([sampletrig(i) sampletrig(i)],y)
%             end
% %             figure,plot(detTrig)
% %             hold on
%             %             plot(resp, 'r')
%             k = 1; rs = []; ws =[];
%             for i=1:length(sampletrig)-1
%                 y = ylim; % current y-axis limits
%                 idx1 = sampletrig(i):sampletrig(i+1);
%                 tmp = difResp(idx1);
%                 tmp2 = detTrig(idx1);
%                 idx = find(tmp > 100);
%                 if length(idx) == 2
%                     rs(k) = idx1(idx(1)); % resposne sample
% %                     plot([rs(k) rs(k)],y);
%                     %                     plot([sampletrig(i) sampletrig(i)],y)
%                     %                     rs(k) - sampletrig(i)
%                     idx2 = find(tmp2 == 8);
%                     if idx1 > 1
%                         ws(k) = idx1(idx2(2)); % waiting sample
% %                         plot([ws(k) ws(k)],y);
% %                         rs(k) - sampletrig(i)
%                     end
%                     k=k+1;
%                 end
%             end
%             %- Reaction time
%             figure,
%             bar((rs - ws)/1e3, 0.4); % (resposne - waiting) / fs
%             set(gca,'Xtick', 1:length(ws),'XtickLabel',1:length(ws));
%             title(['Reaction time, mean:', num2str(mean((rs - ws)/1e3)), ' sec'])
%             xlabel('Trials'), ylabel('time (sec)'),
%             grid on
%             set(gca,'color','none');
%             disp([num2str(mean((rs - ws)/1e3)), ' +- ', num2str(std((rs - ws)/1e3))]);
            
            %%
            
            %%
            disp('preprocessing ...');
            cfg = [];
            cfg.eventid = val1;
            cfg.epochtype = epoch_type;
            cfg.datafile  = datafile;
            cfg.hpfreq = 0.1;
            cfg.lpfreq = 40;
            [f_data] = vy_preprocess(cfg);
            disp('filtering was completed');
%             save(savepath, 'f_data', '-v7.3');
        end
    end
    
    %% Bad trials & channels (automated)
    if flag.preprocessing.artifact == 1
        %
        savepath = fullfile(outd.sub,['a_',subj,'.mat']);
        if exist(savepath, 'file') == 2
            load(savepath)
        else
            cfg                         = [];
            cfg.dataset                 = datafile;
            cfg.trialfun                = 'ft_trialfun_general'; % this is the default
            cfg.trialdef.eventtype      = epoch_type;
            cfg.trialdef.eventvalue     = val1; % the value of the stimulus trigger for fully incongruent (FIC).
            cfg.trialdef.prestim        = 0; % in seconds
            cfg.trialdef.poststim       = 1; % in seconds
            trl = ft_definetrial(cfg);
            
            
            D = trl.trl(:,1);
            [C,idx] = intersect(D-1,sampletrig);
            [C1,idx1] = intersect(sampletrig,D-1);
            Condition1 = Condition(idx1);
            
            cfg = [];
            cfg.trials = idx;
            f_data_spendl = ft_selectdata(cfg,f_data);
            
            %-
            cfg = [];
            cfg.toilim = [-.4,2];
            f_data_spendl = ft_redefinetrial(cfg, f_data_spendl);
            data_in = f_data_spendl;
            cfg = [];
            cfg.pflag = 1; % yes:1, No:2
            cfg.saveflag = 1; % yes:1, No:2
            cfg.savepath = savepath;
            cfg.latency = [data_in.time{1}(1),data_in.time{1}(end)];%[-200,900]./1000;
            cfg.rejectpercentage = .95;
            [r_data,report] = vy_artifactreject(cfg, data_in);
            % disp('Bad data rejection was completed');
        end
        condition_update = Condition1;
        condition_update(report.btrl) = [];
    end
    %%
    if flag.preprocessing.ica == 1
        cfg = [];
        cfg.savepath = outd.sub;
        cfg.savefile = fullfile(outd.sub,['ic_',subj,'.mat']);
        cfg.saveflag = 1;
        cfg.overwrite = 2;
        cfg.lay = lay;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = ic_selection;
        %         cln_data = vy_ica_cleaning_light(cfg, r_data);
        [cln_data,report_ic] = vy_ica_cleaning_light2(cfg, r_data);
        disp('ICA cleaning was completed');
    end
    condition_update(report_ic.btrl) = [];
    
    %%
    savepath = fullfile(outd.sub,['cond_',subj,'.mat']);
    %     condinfo = [];
    %     condinfo.report = report;
    %     %     condinfo.idx = idx;
    %     condinfo.report_ic = report_ic;
    save(savepath, 'condition_update');
        
end

