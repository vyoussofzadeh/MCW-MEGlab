clear; clc, close('all'); warning off

%% Flags
flag.preprocessing.filtering = 1;
flag.preprocessing.artifact = 1;
flag.preprocessing.ica = 1;
flag.notch = 1;
flag.freq = 1;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis
flag.anatomy = 1;  % grand average analysis
flag.sourceanalysis = 1;     % grand average analysis
flag.warping = 1;  % warping to a template, not recommended for clinical reports

%% Initial settings
% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')


cd '/MEG_data/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/MEG_data/epilepsy';
%- Output dir
outdir = '/MEG_data/Vahab/Processed_data';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/MEG_data/Vahab/Github/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
cd('/MEG_data/Vahab/Processed_data/ft_process');
[subjdir] = uigetdir;

%%
disp('1: Definition naming')
disp('2: Picture naming');
task = input('Eneter the task: ');
switch task
    case 1
        %- Auditory definition naming
        tag = 'DFN';
    case 2
        %- Visual picture naming
        tag = 'PN';
end

%%
addpath(allpath.connpath);
addpath(allpath.spm_path);

%%
cd(fullfile(subjdir, tag))

%%
% d = rdir(tag'*22Hz.mat');
d = rdir([subjdir,['/',tag,'/dics/*20Hz.mat']]);

clear datafile
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafile{i} = d(i).name;
end
datafile1 = datafile';
disp(datafile1)
disp('============');

%%
load(datafile1{1});

%%
outputdir = 'results';
if exist(outputdir, 'file') == 0, mkdir(outputdir), end

%%
savenii = fullfile(outputdir,[tag,'_dics.nii']);

template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
cfg = [];
cfg.mask = 'pow';
% cfg.loc = 'min';
cfg.template = template_mri;
cfg.savefile = [];
cfg.volnorm     = 2; % yes: 1
D1 = vy_source_plot(cfg, source_diff_dics);
vy_savenifti(D1,'pow',savenii);


%%
source = ft_read_mri(savenii);
projthresh = 0.65;
s_vol = vy_vol_thresh(source, projthresh, 'anatomy'); % abs

Opt = [];
Opt.savenii = 1; Opt.savefig = 0;
Opt.savename = fullfile(outputdir,['dics_',tag]);
vy_surfce_vis2(s_vol,savenii, Opt);
