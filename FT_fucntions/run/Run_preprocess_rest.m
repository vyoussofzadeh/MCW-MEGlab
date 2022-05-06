
if exist(fullfile(outd.sub,['ic_',subj,'.mat']), 'file') == 2
    load(fullfile(outd.sub,['ic_',subj,'.mat']));
else
    
    if flag.preprocessing.filtering == 1
        %% Filteting, Event reading, Artifact rejecrtion
        savepath = fullfile(outd.sub,['f_',subj,'.mat']);
        if exist(savepath, 'file') == 2
            load(savepath)
        else
            %% Preprocessing
            savepath = fullfile(outd.sub,['f_',subj,'.mat']);
            cfg         = [];
            cfg.channel = {'MEG'};
            cfg.datafile  = datafile;
            f_data        = ft_preprocessing(cfg);
            
            % split dat into 1 s segments
            cfg         = [];
            cfg.length  = 4;
            cfg.overlap = 0;
            dat         = ft_redefinetrial(cfg, f_data);
            dat.time(:) = {(1:size(dat.trial{1},2))/dat.fsample};
            
            % remove the DC-component
            cfg        = [];
            cfg.demean = 'yes';
            cfg.dftfilter = 'yes';
            cfg.hpfilter = 'yes';
            cfg.lpfilter = 'yes';
            cfg.hpfiltord = 3;
            cfg.lpfreq = 70;
            cfg.hpfreq = 0.1;
            cfg.dftfreq = 50;
            cfg.channel = {'megmag', 'meggrad'};
            dat        = ft_preprocessing(cfg, dat);
            
            cfgrsmp = [];
            cfgrsmp.resamplefs  = 500;
            cfgrsmp.detrend     = 'no';
            f_data      = ft_resampledata(cfgrsmp, dat);
            
%             disp('filtering was completed');
%             save(savepath, 'f_data', '-v7.3');
        end
    end
    
    %% Bad trials & channels (manual)
%     clear r_data
%     cfg = [];
%     cfg.metric = 'zvalue';  % use by default zvalue method
%     cfg.layout   = lay;   % this allows for plotting individual trials
%     r_data   = ft_rejectvisual(cfg, f_data);
    
    %% Bad trials & channels (automated)
    if flag.preprocessing.artifact == 1
        
        savepath = fullfile(outd.sub,['a_',subj,'.mat']);
        cfg = [];
        cfg.pflag = 2; % yes:1, No:2
        cfg.saveflag = 2; % yes:1, No:2
        cfg.savepath = savepath;
        cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)];%[-200,900]./1000;
        cfg.rejectpercentage = .95;
        [r_data,report] = vy_artifactreject(cfg, f_data);
    end
    
    
    %% ICA cleaning
    if flag.preprocessing.ica == 1
        cfg = [];
        cfg.savepath = outd.sub;
        cfg.savefile = fullfile(outd.sub,['rest_ic_',subj,'.mat']);
        cfg.saveflag = 1;
        cfg.overwrite = 2;
        cfg.lay = lay;
        cfg.run = run;
        cfg.n   = 20;
        cfg.subj = subj;
        cfg.allpath = allpath;
        cfg.select = 1;
        %         cln_data = vy_ica_cleaning_light(cfg, r_data);
%         cln_data = vy_ica_cleaning_light3(cfg, r_data);
        cln_data = vy_ica_cleaning_light6(cfg, r_data);
        disp('ICA cleaning was completed');
    end
    
    %%
    %     cfgrsmp = [];
    %     cfgrsmp.resamplefs  = 500;
    %     cfgrsmp.detrend     = 'no';
    %     f_data      = ft_resampledata(cfgrsmp, f_data);
    %
    %     %% ICA cleaning
    %     cfg = [];
    %     cfg.savepath = fullfile(outd.sub,['ic_',subj,'.mat']);
    %     cfg.saveflag = 1;
    %     cfg.lay = lay;
    %     cfg.n   = 20;
    %     cln_data = vy_ica_cleaning_light(cfg, f_data);
    %     disp('ICA cleaning was completed');
    %
    %     %%
    %     % muscle
    %     cfg            = [];
    %     %     cfg.trl        = trl;
    %     %     cfg.datafile   = 'ArtifactMEG.ds';
    %     %     cfg.headerfile = 'ArtifactMEG.ds';
    %     %     cfg.continuous = 'yes';
    %
    %     % channel selection, cutoff and padding
    %     cfg.artfctdef.zvalue.channel      = 'meg*';
    %     cfg.artfctdef.zvalue.cutoff       = 4;
    %     cfg.artfctdef.zvalue.trlpadding   = 0;
    %     cfg.artfctdef.zvalue.fltpadding   = 0;
    %     cfg.artfctdef.zvalue.artpadding   = 0.1;
    %
    %     % algorithmic parameters
    %     cfg.artfctdef.zvalue.bpfilter     = 'yes';
    %     cfg.artfctdef.zvalue.bpfreq       = [110 140];
    %     cfg.artfctdef.zvalue.bpfiltord    = 9;
    %     cfg.artfctdef.zvalue.bpfilttype   = 'but';
    %     cfg.artfctdef.zvalue.hilbert      = 'yes';
    %     cfg.artfctdef.zvalue.boxcar       = 0.2;
    %
    %     % make the process interactive
    %     cfg.artfctdef.zvalue.interactive = 'yes';
    %
    %     [cfg, artifact_muscle] = ft_artifact_zvalue(cfg, cln_data);
    %
    %     %% Converting to fif file
    %     isepch = ft_datatype(cln_data, 'raw');
    %
    %     hdr = ft_read_header(datafile);
    %     test_var = cln_data;
    %     test_var.hdr = hdr;
    %     fieldtrip2fiff('test.fif', test_var)
    
    %%
    
    %     %%
    %     cfg                         = [];
    %     cfg.dataset                 = datafile;
    %     cfg.channel = {'MEG'};
    %     f_data = ft_preprocessing(cfg);
    %
    %     %%
    %     cfgrsmp = [];
    %     cfgrsmp.resamplefs  = 500;
    %     cfgrsmp.detrend     = 'no';
    %     f_data      = ft_resampledata(cfgrsmp, f_data);
    %
    %     %%
    %     cfg            = [];
    %     cfg.method     = 'runica';
    %     cfg.numcomponent = 20;       % specify the component(s) that should be plotted
    %     comp           = ft_componentanalysis(cfg, f_data);
    %
    %     %%
    %     cfg           = [];
    %     cfg.component = 1:20;       % specify the component(s) that should be plotted
    %     cfg.layout    = lay;
    %     cfg.comment   = 'no';
    %     ft_topoplotIC(cfg, comp)
    %     colormap(brewermap(256, '*RdYlBu'));
    %     %%
    %     cfg = [];
    %     cfg.viewmode = 'component';
    %     cfg.layout = lay;
    %     ft_databrowser(cfg, comp);
    %     colormap(brewermap(256, '*RdYlBu'));
    %     % set(gcf, 'Position', [600   600   700   500]);
    %
    %     %%
    %     cfg         = [];
    %     cfg.length  = 2;
    %     cfg.overlap = 0;
    %     dat         = ft_redefinetrial(cfg, f_data);
    
    %% Bad channels - trials
    %     savepath = fullfile(outd.sub,['a_',subj,'.mat']);
    %     cfg = [];
    %     cfg.pflag = 2; % yes:1, No:2
    %     cfg.saveflag = 1; % yes:1, No:2
    %     cfg.savepath = savepath;
    %     cfg.latency = [dat.time{1}(1),dat.time{1}(end)];%[-200,900]./1000;
    %     [r_data,report] = vy_artifactreject_rest(cfg, dat);
    %     % disp('Bad data rejection was completed');
    
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
    %     cfg.continuous = 'yes';
    %     %     cfg.trials     = report.btrl;
    %     cfg.channel   = report.bchan;
    %     ft_databrowser(cfg,f_data);
    %
    %     %% ICA cleaning
    %     cfg = [];
    %     cfg.savepath = fullfile(outd.sub,['ic_',subj,'.mat']);
    %     cfg.saveflag = 1;
    %     cfg.lay = lay;
    %     cfg.n   = 20;
    %     cln_data = vy_ica_cleaning_light(cfg, r_data);
    %     disp('ICA cleaning was completed');
    
    %%
    %     tmpoptions   = {'doplot', 'no', 'trial', [] , 'plottype', 'summary'}; % perform frequency analysis
    %     comp_freq = hcp_ICA_freq(cln_data.comp, tmpoptions);
    
    %%
    
    
    
    
    
end