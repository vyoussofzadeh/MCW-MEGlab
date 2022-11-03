%% BCI dataset, Ulster University & Medical College of Wisconsin

% Script: Preprocess (band pass filtering, ICA, reject trials)
% Project: BCI_nonstationarity
% Writtern by: Vahab Youssof Zadeh
% Update: 11/03/2022

clear; clc, close('all'); warning off

%% Paths
addpath('./data')
addpath('./run')
Run_setpath

%% Listing data
cd(indir)

clear d_all sub session
d = rdir(fullfile(indir,'/**/*meg.mat'));
for i=1:length(d)
    tkz = tokenize(d(i).name,'/');
    D.df{i} = tkz{end};
    
    tkz1 = tokenize(D.df{i},'_');
    tkz2 = tokenize(tkz1{1},'-');
    
    session{i} = tkz1{2}(end);
    sub{i} = tkz2{2};
    
    D.ss{i} = [sub{i}, '_', session{i}];
    D.sub{i} = sub{i};
    D.session{i} = session{i};
    D.datafile{i} = d(i).name;
end
disp(D.ss') % Data subject/Session
disp(D.sub') % Data subject/Session

%%
D_sel=  [];
for i=1:length(D.ss)
    D_sel{i} = [num2str(i), ':', D.ss{i}];
end
disp(D_sel') % Data subject/Session
data_sel = input('sel data: ');

%% Read and Preprocess
for df = 1:length(data_sel)
    
    d_sel = data_sel(df);
    subjdir = D.ss{d_sel};
    datafile = D.datafile{d_sel};
    
    %% 4D layout
    cfg = [];
    cfg.layout = 'neuromag306mag.lay';
    lay = ft_prepare_layout(cfg);
    % ft_layoutplot(cfg);
    disp('============');
    
    %%
    disp(datafile)
    disp(['subj-session:',D.ss{d_sel}])
    %     disp(['Tag:',data_tag])
    disp('============');
    
    subj = [D.sub{d_sel}, '_', D.session{d_sel}];
    
    %%
    meg_data = load(datafile);
    trl = meg_data.dataMAT.trialinfo;
    meg_data = meg_data.dataMAT;
    
    %%
    cfg = [];
    cfg.toilim = [0,5];
    ep_data = ft_redefinetrial(cfg, meg_data);
    
    %%
    outd.sub = fullfile(outdir, D.sub{d_sel}, D.session{d_sel});
    if exist(outd.sub, 'file') == 0
        mkdir(outd.sub);   %create the directory
    end
    cd(outd.sub)
    disp(['outputdir:',outd.sub])
    disp('============');
    
    savefile = fullfile(outd.sub,['ic_',D.sub{d_sel}, '_', D.session{d_sel}, '.mat']);
    if exist(savefile, 'file') == 2
        disp('already saved!')
    else
        
        %- Preprocessing, bandpass filtering
        cfg = [];
        cfg.hpfilter = 'yes';
        cfg.lpfilter = 'yes';
        cfg.dftfilter = 'yes';
        cfg.hpfiltord = 3;
        cfg.hpfreq = 1;
        cfg.lpfreq = 40;
        cfg.channel = {'MEG'};
        f_data = ft_preprocessing(cfg, ep_data);
        
        %- Preprocessing, rejecting bad trials
        cfg = [];
        cfg.pflag = 2; % yes:1, No:2
        cfg.saveflag = 0; % yes:1, No:2
        cfg.savepath = [];
        cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)];
        cfg.rejectpercentage = .95;
        cfg.method =  'manual'; %'manual'; % 'auto';
        [r_data,~] = do_rejectdata(cfg, f_data);
        
        %- ICA preprocessing
        cfg = [];
        cfg.lay = lay;
        cfg.subj = subj;
        cfg.n = 20;
        cfg.allpath = allpath;
        cfg.savefig = 1;
        data_ica = do_ica_bci(cfg, r_data);
%         data_ica.grad = []; data_ica.cfg = [];
        
        close all
        %- Save preprocessed data
        disp('saving data ...')
        
        cd(outd.sub)
        save(savefile, 'data_ica', '-v7.3');
        
        addpath(genpath(allpath.script_path));
        addpath(allpath.ft_path); ft_defaults
        
    end
end

