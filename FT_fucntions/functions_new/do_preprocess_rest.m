function cln_data = do_preprocess_rest(cfg_main, datafile)

subdir = cfg_main.datalog.subdir;
subj = cfg_main.datalog.subj;
flag = cfg_main.flag;
lay = cfg_main.lay;
all_path = cfg_main.all_path;
run = cfg_main.datalog.run;

%%
savefile = fullfile(subdir,['ic_',subj,'_run_', run, '.mat']);
if exist(savefile, 'file') == 2
    load(savefile);
else
    
    if flag.preproces.filtering == 1
        %% Filteting, Event reading, Artifact rejecrtion
%         savepath = fullfile(subdir,['f_',subj,'.mat']);
%         if exist(savepath, 'file') == 2
%             load(savepath)
%         else
            %% Preprocessing
%             savepath = fullfile(subdir,['f_',subj,'.mat']);
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
            cln_data      = ft_resampledata(cfgrsmp, dat);
            
            %             disp('filtering was completed');
            %             save(savepath, 'f_data', '-v7.3');
%         end
    end
    
    %% Bad trials & channels (manual)
    %     clear r_data
    %     cfg = [];
    %     cfg.metric = 'zvalue';  % use by default zvalue method
    %     cfg.layout   = lay;   % this allows for plotting individual trials
    %     r_data   = ft_rejectvisual(cfg, f_data);
    
    %% Bad trials & channels (automated)
    if flag.preproces.artifact == 1
        cfg = [];
        cfg.pflag = 2; % yes:1, No:2
        cfg.saveflag = 0; % yes:1, No:2
        cfg.savepath = [];
        cfg.latency = [cln_data.time{1}(1),cln_data.time{1}(end)];
        cfg.rejectpercentage = .95;
        cfg.method =  'auto'; %'manual'; % 'auto';
        [cln_data,~] = do_rejectdata(cfg, cln_data);
    end
    
    
    %% ICA cleaning
    if flag.preproces.ica == 1
        
        cfg = [];
        cfg.lay = lay;
        cfg.subj = [subj,'_', run];
        cfg.n = 20;
        cfg.allpath = all_path;
        cfg.savefig = 1;
        cln_data = do_ica(cfg, cln_data);
        disp('ICA cleaning was completed');
    end
    
    close all
    %- Save preprocessed data
    disp('saving data ...')
    
    cd(subdir)
    save(savefile, 'cln_data', '-v7.3');
    
    addpath(genpath(all_path.script_path));
    addpath(all_path.ft_path); 
    ft_defaults
    
end
end