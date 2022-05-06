%%
clear; clc, close('all'); warning off

%% Flags
% flag.preprocessing.filtering = 1;
% flag.preprocessing.artifact = 1;
% flag.preprocessing.ica = 1;
% flag.notch = 1;
% flag.freq = 1;     % TFR & FFT
% flag.time = 1;     % Time-locked & Cov estimation
% flag.gave = 0;     % grand average analysis
% flag.anatomy = 1;     % grand average analysis
% flag.sourceanalysis = 1;     % grand average analysis

%% Initial settings
set(0,'DefaultFigureWindowStyle','docked')

cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));
rmpath('./Failedattemps');


%- Input dir
indir = '/data/MEG/Clinical/MEG';
%- Output dir
outdir = '/data/MEG/Clinical';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
msk = 'pow';
projthresh = 0.6;

addpath(allpath.connpath);
addpath(allpath.spm_path);
close all

%%
cd ('/data/MEG/Clinical/ft_process/');
[subdir] = uigetdir;
% subdir = '/data/MEG/Clinical/ft_process/17/alby_tracy';
% subdir = '/data/MEG/Clinical/ft_process/19/bednar_peggy';

%%
disp('1: Definition naming')
disp('2: Picture naming');
disp('3: both')
task = input('Eneter the task: ');

%%
switch task
    case 1
        tsk1 = {'DFN'};
    case 2
        tsk1 = {'PN'};
    case 3
        tsk1 = {'DFN', 'PN'};
end

for i=1:length(tsk1)
    
    tsk = tsk1{i};
    %%
    datadir = fullfile(subdir,tsk, 'dics');
    Index = strfind(datadir, '/'); sub  = datadir(Index(6)+1:Index(7)-1);
    
    file = rdir([datadir '/*.mat']);
    files = {file.name}';
    files_date = {file.date}';
    [dn,idx] = sort(datenum(files_date, 'dd-mm-yyyy hh:MM:ss'), 1, 'descend');
    
    
    filesaved = files{idx(1)};
    
    % filesaved = 'dics_alby_tracy_17Hz.mat'; sub = 'bednar_peggy';
    
    cd(subdir)
    load(filesaved)
    
    %%
    % [subjdir] = uigetdir;
    % subdir = '/data/MEG/Clinical/ft_process/17/alby_tracy/PN/dics';
    % cd(subdir)
    % filesaved = 'dics_alby_tracy_16Hz.mat';
    % sub = 'alby_tracy';
    % tsk = 'PN';
    % load(fullfile(subdir,filesaved))
    %%
    
    cfg = [];
    cfg.mask = 'pow';
    cfg.loc = 'min';
    cfg.template = template_mri;
    cfg.savefile = [];
    cfg.volnorm     = 2; % yes: 1
    source_dics = vy_source_plot(cfg, source_diff_dics);
    
    
    projthresh = 0.6;
    source_dics1 = vy_vol_thresh(source_dics,projthresh,'pow'); % abs
    source_dics1.pow = -(source_dics1.pow./max(source_dics1.pow(:)));
    
    set(gcf,'Name',sub) %select the name you want
    savenname = tsk;
    vy_savenifti(source_dics1,'pow',[savenname,'.nii']);
    %     print(['./', tsk,'/',names{i}],'-depsc');
    %     print(savenname,'-dpng');
    
    %-Surf-vis
    opt.run = [];
    opt.tsk = tsk;
    opt.subj = sub;
    opt.savedir = [];
    opt.savenii = 2;
    %     opt.plot = '-mosaic';
    opt.plot = '-row';
    vy_surfce_vis([],[savenname,'.nii'], opt);
    
end
