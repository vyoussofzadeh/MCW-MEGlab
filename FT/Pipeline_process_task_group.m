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
cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/data/MEG/Clinical/MEG';
%- Output dir
outdir = '/data/MEG/Clinical';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
disp('1: Definition naming')
disp('2: Picture naming');
task = input('Eneter the task: ');

switch task
    case 1
        % - Auditory definition naming
        tsk = 'DFN'; tag1 = 'dfn';
    case 2
        % - Visual picture naming
        tsk = 'PN';
end

DestDirectory = fullfile(outdir);

d = rdir([DestDirectory,['/**/',tsk,'/dics/*22Hz.mat']]);
% d = rdir([DestDirectory,['/**/19/**/',tsk,'/dics/*22Hz.mat']]);
% d = rdir([DestDirectory,['/**/18/**/',tsk,'/dics/*22Hz.mat']]);
% d = rdir([DestDirectory,['/**/14_5/**/',tsk,'/dics/*22Hz.mat']]);


% d = rdir([DestDirectory,['/**/',tsk,'/dics/*.mat']]);
% d = rdir([DestDirectory,['/19/**/',tsk,'/dics/*.mat']]);
% d = rdir([DestDirectory,['/18/**/',tsk,'/dics/*.mat']]);
% d = rdir([DestDirectory,['/17/**/',tsk,'/dics/*.mat']]);
% d = rdir([DestDirectory,['/16/**/',tsk,'/dics/*.mat']]);
% d = rdir([DestDirectory,['/14_5/**/',tsk,'/dics/*.mat']]);
% d = rdir([DestDirectory,['/13/**/',tsk,'/dics/*.mat']]);
% d = rdir([DestDirectory,['/12/**/',tsk,'/dics/*.mat']]);
% d = rdir([DestDirectory,['/11/**/',tsk,'/dics/*.mat']]);

files = {d.name}'; 
files_sel = files(1:end,1);

%%
outputdir = fullfile(DestDirectory,'group_dics');
if exist(outputdir, 'file') == 0
    mkdir(outputdir);   %create the directory
end
cd(outputdir)

%% Group source analysis (averaging)
Run_group_analysis

%% Group visualisation
Run_group_vis

%% Task comparison (PN vs DFN)
Run_taskcompare


