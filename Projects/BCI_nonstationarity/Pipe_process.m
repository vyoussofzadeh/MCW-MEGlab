%% BCI dataset, Ulster University & Medical College of Wisconsin

% Script: Process (extract virtual sensrors at voxel and roi AAL atlas levels)
% Project: BCI_nonstationarity
% Writtern by: Vahab Youssof Zadeh
% Update: 11/03/2022

clear; clc, close('all'); warning off

%% Paths
addpath('./run')
Run_setpath
addpath('./data')
addpath('./run')
addpath(genpath('./functions'))

%% Loading up data
cd(outdir)

clear d_all sub session
d = rdir(fullfile(outdir,'/**/ic_*.mat'));
for i=1:length(d)
    tkz = tokenize(d(i).name,'/');
    D.df{i} = tkz{end};
    
    tkz1 = tokenize(D.df{i},'_');
    tkz2 = tokenize(tkz1{end},'.mat');
    
    sub{i} = tkz1{2};
    session{i} = tkz2{1};
    
    D.ss{i} = [sub{i}, '_', session{i}];
%     disp(D.ss{i})
    D.sub{i} = sub{i};
    D.session{i} = session{i};
    D.datafile{i} = d(i).name;
end
disp(D.ss') % Data subject/Session

%%
D_sel=  [];
for i=1:length(D.ss)
    D_sel{i} = [num2str(i), ':', D.ss{i}];
end
disp(D_sel') % Data subject/Session
data_sel = input('sel data: ');

%% anatomy (template) 
load('anat')

%% Preprocess
for df = 1:length(data_sel)
    
    d_sel = data_sel(df);
    subj = D.ss{d_sel};
    disp(['Analyzing: ', subj])
    
    idx = strfind(D.datafile{d_sel},'/');
    D.datafile{d_sel}(1:idx(end)-1);
    
    savefile = fullfile(D.datafile{d_sel}(1:idx(end)-1),['vs_', D.datafile{d_sel}(idx(end)+4:end-4),'.mat']); % output dir
    if exist(savefile, 'file') == 2
        disp('already saved!')
    else
        load(D.datafile{d_sel})
        
        %- epoching
        cfg = [];
        cfg.toilim = [0,5];
        ep_data = ft_redefinetrial(cfg, data_ica);
        
        %- Source analysis, time-domain beamformer (LCMV)
        cfg = [];
        cfg.individual_grid = individual_grid;
        cfg.vol = individual_headmodel;
        source_active = do_sourceanalysis(cfg, data_ica);
        
        cfg = [];
        cfg.individual_grid = individual_grid;
        cfg.atlas = atlas;
        [vs, vs_roi] = do_extractvirtualsensor(cfg, source_active); %- Extract virtual sensors
        
        trialinfo = data_ica.trialinfo;
        save(savefile, 'vs_roi', 'vs', 'trialinfo','-v7.3');
        
    end
end

