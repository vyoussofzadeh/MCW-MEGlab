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

%% check if speech is available
hdr = ft_read_header(datafile);
Index = strfind(hdr.label,{'MISC001'});
Index = find(not(cellfun('isempty',Index)));
disp(Index),


if ~isempty(Index)
    
    %% Preprocesssing
    clear cln_data
    ic_selection = 1; % 1: manual, 2: automated
    switch task
        
        case {1,2}
            Run_preprocess
            %- Speech
            if flag.speechanalysis ==1
                Run_speech
            end
    end
    disp('no speech data available')
    
end




