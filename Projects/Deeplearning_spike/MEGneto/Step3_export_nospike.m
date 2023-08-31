%% The Spike Detection MEG pipline

% Spike detection MEG pipline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 08/31/2023

clear; clc, close('all'); warning off

%% Flags
flag.preprocessing.filtering = 1;
flag.preprocessing.artifact = 1;
flag.preprocessing.ica = 1;
flag.notch = 1;
flag.freq = 1;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis
flag.anatomy = 1;     % grand average analysis
flag.sourceanalysis = 1;     % grand average analysis
flag.speechanalysis = 1;     % speech analysis
flag.analysis = 1;

%% Datalog (subject details)
Datalog = [];

%% Initial settings
cd '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/MEG_data/Research_studies/Epil_annotated_data/annotated_info';
%- Output dir
outdir = '/MEG_data/Research_studies/Epil_annotated_data/annotated_data_nospike';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/MEG_data/LAB_MEMBERS/Vahab/Github/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
cd(indir)
d = rdir([indir,'/*.mat']);

%%
clear subj run sub_run
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    tkz = tokenize(name,'_');
    subj{i} = [tkz{1}, '_', tkz{2}];
    sub_run{i,:} = [subj{i}, '_', tkz{3}];
end
[sub_run_unq,IA,IC] = unique(sub_run);
disp(sub_run_unq);

[sub_unq,IA_s,IC_s] = unique(subj);
disp(sub_unq);

%% Read, preprocess, and save data
if exist(outdir, 'file') == 0, mkdir(outdir), end

for i=1:length(d)
    disp([num2str(i),'/',num2str(length(d))])
    [pathstr, name] = fileparts(d(i).name);
    load(d(i).name);
    if ~exist([fullfile(outdir,name), '.mat'], 'file') && exist(Anot.filename, 'file')
        %     datafile = Anot.filename;
        disp('preprocessing ...')
        cfg = [];
        cfg.dataset = Anot.filename;
        cfg.channel = {'megmag', 'meggrad', 'eeg'};
        raw_data = ft_preprocessing(cfg);
        
        cfg = [];
        cfg.resamplefs = 500;
        dsample_data = ft_resampledata(cfg, raw_data);
        
        anot_data_all = []; k=1;
        for j=1:size(Anot.T_tint,1)
            cfg = [];
            if  Anot.T_tint(j,2) -  Anot.T_tint(j,1) > 0.2 && Anot.T_tint(j,1) > 5
                cfg.toilim = [Anot.T_tint(j,1)-2, Anot.T_tint(j,1)-1.6]./(raw_data.fsample/dsample_data.fsample); % no spike
                anot_data = ft_selectdata(cfg,dsample_data);
                anot_data_all{k} = anot_data;
                k=k+1;
            end
        end
        tkz = tokenize(name,'_');
        save(fullfile(outdir,name), 'anot_data_all', 'Anot')
    end
end
