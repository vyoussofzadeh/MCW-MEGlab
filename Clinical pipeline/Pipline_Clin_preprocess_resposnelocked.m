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
disp('1: Definition naming')
disp('2: Picture naming');
disp('3: Motor');
task = input('Eneter the task: ');
switch task
    case 1
        %- Auditory definition naming
        tag = 'DFN';
    case 2
        %- Visual picture naming
        tag = 'PN';
    case 3
        %- Motor task
        tag = 'motor';
end

%%
cd(indir)
[subjdir] = uigetdir;

%%
d = rdir([subjdir,['/**/','sss','/*',tag,'*/*raw_tsss.fif']]);

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
   datasel = input('choose which data to analyze, row number:'); 
else
    datasel = 1;
end
disp([subj, ' and,'])
disp([datafile1{datasel}, 'was selected for the analysis ...'])
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
Date  = datafile(Index(4)+1:Index(5)-1);
disp('============');
disp(datafile)
disp(['subj:',subj])
disp(['Date:',Date])
disp('============');

%%
%-elec/grad
sens = ft_read_sens(datafile);
sens = ft_convert_units(sens,'mm');

%%
outd.sub = fullfile(outdir,'ft_process',subj, tag);
if exist(outd.sub, 'file') == 0
    mkdir(outd.sub);   %create the directory
end
cd(outd.sub)
disp(['outputdir:',outd.sub])
disp('============');

%% Preprocesssing
clear cln_data
ic_selection = 1; % 1: manual, 2: automated
switch task
    
    case {1,2}
        Run_preprocess
        
        %- Notch filtering: 30Hz
        if flag.notch ==1, Run_notch, end %- Notch filtering of 30Hz
        
        savetag1 = [tag,'_',subj,'_BS_channels_reslocked'];
        savetag2 = [tag,'_',subj,'_IC_data_reslocked'];
        %- BrainStorm export preprocessed ft-IC
        
        %- Speech
        if flag.speechanalysis ==1
            Run_speech
        end
        
        %%
        % export to brainstorm
        Run_bs
        
    case 3
        Run_preprocess_motor
        cln_data = clnl_data; % left data condition
        if flag.notch ==1, Run_notch, end %- Notch filtering of 30Hz
        savetag1 = [tag,'_',subj,'_left_BS_channels'];
        savetag2 = [tag,'_',subj,'_left_IC_data'];
        Run_bs
        
        cln_data = clnr_data; % right data condition
        if flag.notch ==1, Run_notch, end
        savetag1 = [tag,'_',subj,'_right_BS_channels'];
        savetag2 = [tag,'_',subj,'_right_IC_data'];
        Run_bs
end

%% Selecting freq of interest
fmax = 40;
datain = cln_data;
if flag.freq == 1
    stag = 'tsk_baseline';
    Run_freq
    disp(['time_of_interest:',num2str(time_of_interest),'sec']);
    disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
    L = 0.3;
end
pause(10)

%%
close all
set(0,'DefaultFigureWindowStyle','normal')
addpath('/usr/local/MATLAB_Tools/brainstorm3')
%brainstorm

%%
cd (bssavedir)
clear

disp('completed, data are ready to import into BS!')




