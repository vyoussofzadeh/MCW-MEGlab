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
flag.speechanalysis = 2;     % speech analysis

%% Initial settings
set(0,'DefaultFigureWindowStyle','docked')

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
tag = 'spont';

%%
cd(indir)
[subjdir] = uigetdir;
cd(subjdir)

%%
d = rdir([subjdir,['/**/','sss','/*',tag,'*/*raw_tsss.fif']]); stag = 'tsss_';
% d = rdir([subjdir,['/**/','sss','/*',tag,'*/*t_sss_ecgClean_raw.fif']]); stag = 'sss_ecgClean_raw';

%%
clear subj datafolder datafile datafile1
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    subj = datafile{i}(Index(3)+1:Index(4)-1);
end
datafile1 = datafile';
disp(datafile1)
if length(datafile1) > 1
    datasel = input('choose data to analyze, eg, 1,2:');
else
    datasel = 1;
end
disp([subj, ' and,'])
disp(datafile1{datasel})
disp('was selected for the analysis.')
disp('============');

%%
epoch_type = 'STI101';

%% 4D layout
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);
disp('============');

%%
close all
datafile = datafile1{datasel}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
Index = strfind(datafile, '/');
% Date  = datafile(Index(4)+1:Index(5)-1);
Run   = datafile(Index(end)+4:Index(end)+5);
disp('============');
disp(datafile)
disp(['subj:',subj])
% disp(['Date:',Date])
disp(['Run:',Run])
disp('============');

%%
%-elec/grad
sens = ft_read_sens(datafile);
sens = ft_convert_units(sens,'mm');

%%
outd.sub = fullfile(outdir,'ft_process',subj, tag, Run);
if exist(outd.sub, 'file') == 0
    mkdir(outd.sub);   %create the directory
end
cd(outd.sub)
disp(['outputdir:',outd.sub])
disp('============');

%% Preprocesssing
clear cln_data
ic_selection = 1; % 1: manual, 2: automated
Run_preprocess_rest

%- Notch filtering: 30Hz
% if flag.notch ==1, Run_notch, end %- Notch filtering of 30Hz
savetag1 = ['Run_',Run,'_', subj,'_', stag, '_BS_channels'];
savetag2 = ['Run_',Run,'_', subj,'_', stag, 'IC_MEG'];
%- BrainStorm export preprocessed ft-IC

%%
% export to brainstorm
Run_bs

%%
close all
set(0,'DefaultFigureWindowStyle','normal')
addpath('/usr/local/MATLAB_Tools/brainstorm3')
%brainstorm

%%
cd (bssavedir)
clear

disp('completed, data are ready to import into BS!')

%%


