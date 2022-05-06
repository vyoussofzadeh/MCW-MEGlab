clear; clc, close('all'); warning off

%% Flags
flag.preprocessing.filtering = 1;
flag.preprocessing.artifact = 1;
flag.preprocessing.ica = 1;
flag.notch = 1;
flag.freq = 1;         % TFR & FFT
flag.time = 1;         % Time-locked & Cov estimation
flag.gave = 0;         % grand average analysis
flag.sourceanalysis = 1;     % grand average analysis
flag.warping = 1;      % warping to a template, not recommended for clinical reports
flag.speech = 1;
% flag.anatomy = 1;  % anatomy check
flag.anatomycheck = 1; % anatomy check

%% Initial settings
% set(0,'DefaultFigureWindowStyle','docked')
% set(0,'DefaultFigureWindowStyle','normal')

cd('/MEG_data/Vahab/Github/MCW-MEGlab/FT');
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
        %- Visual picture naming
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
sens = ft_read_sens(datafile,'senstype','meg');
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
Run_preprocess

%% Speech
if flag.speech ==1
    if ~exist('./speech', 'dir'), mkdir('speech'), end
    Run_speech
    % Run_speech_2
end

%% Notch filtering: 30Hz
if flag.notch ==1
    Run_notch
end

%% Selecting freq of interest
fmax = 40;
datain = cln_data;
if flag.freq == 1
    stag = 'tsk_baseline';
    Run_freq
    disp(['time_of_interest:',num2str(time_of_interest),'sec']);
    disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
end

%% Response-lock analysis

switch task
    case 1
        Run_timelock_resp
    case 2
        Run_timelock_resp_PN_sel
end

%% Grand averaging
if flag.gave == 1
    Run_grandmean
end

%% Source analysis
if flag.sourceanalysis == 1
    
    outputmridir = fullfile(outdir,'ft_process', subj,'anat'); % output dir
    if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
    %     close all
    
    disp('1: Surface-based')
    disp('2: Volumetric');
    analysis = input('Eneter the analysis: ');
    %     analysis = 2;
    
    switch analysis
        case 1
            mtag = 'source_surf';
            disp('1: SPM source analysis (surface + BF)');
            disp('2: Export ft to Brainstorm');
            disp('3: Source-based bf');
            method = input('Method: ');
        case 2
            % end
            disp('1: LCMV')
            disp('2: Network+Connectvity');
            disp('3: DICS Source');
            method = input('Method: ');
            %             method = 3;
            
            %-
%             disp('1: indiv grid')
%             disp('2: indiv grid, warpped with template')
%             choose_grid = input(': ');
            choose_grid = 2;
            
            %-
            disp('1: Low-res grid')
            disp('2: High-res grid')
            meshgridres = input('Mesh grid: ');
            %             meshgridres = 1;
            
            disp('1: Response-locked');
            disp('2: Time-locked');
            rl = input(':');
            
            disp('1: anatomy check');
            disp('2: no');
            ac = input(':');
            if ac==1
                flag.anatomycheck = 1; % anatomy check
            else
                flag.anatomycheck = 2; % anatomy check
            end            
    end
    
    switch analysis
        case 1
            % Surface-based analysis
            Run_surfacebased
        case 2
            % Volumetric-based analysis
            Run_volumetric
    end
    disp('============');
end
%-
disp([datafile,' ,was completed'])
disp(['output as,', outd.sub])
