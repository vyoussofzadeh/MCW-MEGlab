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
cd('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/FT');
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/MEG_data/epilepsy';
%- Output dir
outdir = '/MEG_data/LAB_MEMBERS/Vahab/Processed_data';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/MEG_data/Vahab/Github/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
disp('1: Definition naming')
disp('2: Picture naming');
disp('3: Motor');
disp('4: Semantic decision');
task = input('Eneter the task: ');
switch task
    case 1
        %- Auditory Definition Naming
        tag = 'DFN';
    case 2
        %- Picture Naming
        tag = 'PN';
    case 3
        %- Motor task
        tg = 'motor';
    case 4
        %- Semantic decision
        tag = 'SD';
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
disp([subj, ' data was selected for the analysis'])
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
datafile = datafile1{1}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
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
Run_preprocess

%% Notch filtering: 30Hz
% if flag.notch ==1
%     Run_notch
% end

%% Selecting freq of interest
datain = cln_data;
fmax = 40;
if flag.freq == 1
    stag = 'tsk_baseline';
    Run_freq
    disp(['time_of_interest:',num2str(time_of_interest),'sec']);
    disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
    L = 0.3;
end

%% Selecting time of interest
% if flag.time == 1
%     %--setting baseline interval
%     toi(1,:) = [-0.3,0];
%
%     %-- setting the post-stim interval
%     disp(['suggested time interval:',num2str(time_of_interest), '+-', num2str(L),' Sec']);
%     %     disp('Yes: 1, No: 2');
%     %     tfa = input('Is it OK to proceed?');
%     %     disp('the following time was selected');
%     %        if (time_of_interest-L) > 0.4 && (time_of_interest+L) < 2.3
%     if (time_of_interest-L) >= 0.1 && (time_of_interest+L) < 2.3
%         toi(2,:) = round([time_of_interest-L, time_of_interest+L],1);
%     else
%         switch task
%             case 1
%                 %                         toi = [-0.3,0;1.1,1.7]; % Best of DN
%                 toi(2,:) = [0.8,1.5]; % Best of DN
%             case 2
%                 toi(2,:) = [0.4,1.2]; % Best for PN, left IFG
%                 toi(2,:) = [0.4,1.6]; % Best for PN, left IFG
%                 toi(2,:) = [0.7,1.6]; % Best for PN, left IFG
%                 toi(2,:) = [1,1.6]; % Best for PN, left IFG
%             case 3
%                 toi(1,:) = [-0.25,0.25];
%                 toi(2,:) = [0.75,1.25]; % Best for PN, left IFG
%         end
%     end
%     disp(['[',num2str(toi(1,1)),',',num2str(toi(1,2)),'] sec interval was selected as bsl']);
%     disp(['[',num2str(toi(2,1)),',',num2str(toi(2,2)),'] sec interval was selected as pst']);
%     Run_time
% end

%%
toi = [-0.3,0;0.4,1.6];
% toi = [-0.3,0;0.4, 0.8];
ep_data = vy_epoch(datain, toi);

cfg = [];
ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);

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
            disp('1: Low-res grid')
            disp('2: High-res grid')
            choose_grid = input('Mesh grid: ');
            meshgridres = 1;
            rl = 2;
    end
    
    switch analysis
        case 1
            savetag1 = [tag,'_',subj,'_BS_channels'];
            savetag2 = [tag,'_',subj,'_IC_data'];
            % Surface-based analysis
            Run_surfacebased
        case 2
            % Volumetric-based analysis
            flag.anatomycheck = 1;
            Run_volumetric
    end
    disp('============');
end
%-
disp([datafile,' ,was completed'])
disp(['output as,', outd.sub])

%%

