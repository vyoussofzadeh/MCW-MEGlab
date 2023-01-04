%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: Process (extract virtual sensrors at voxel and roi AAL atlas levels)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 11/16/2022

clear; clc, close('all'); warning off

%% Paths
addpath('./run')
Run_setpath
addpath('./data')
addpath('./run')
addpath(genpath('./functions'))

%% Loading up data
cd(outdir)

%- Semantic decision
tag = 'SD';

clc
if exist(['datalog_',tag,'.mat'], 'file') == 2
    load(['datalog_',tag,'.mat'])
else
    clear datafolder datafile subj_all
    datafile1 = [];
    d = rdir([indir,['/**/tSSS/*',tag,'*_raw.fif']]);
    for i=1:length(d)
        [pathstr, name] = fileparts(d(i).name);
        datafolder{i} = pathstr;
        datafile{i} = d(i).name;
        Index = strfind(datafile{i}, '/');
        subj_all{i} = [num2str(i), ': ', datafile{i}(Index(6)+1:Index(7)-1)];
    end
    datafile1 = vertcat(datafile1,datafile);
    datafile1 = datafile1';
    save(['datalog_',tag,'.mat'],'datafile1','subj_all')
end
disp(datafile1)
disp('============');

%%
disp('1: choose specific subject');
disp('2: do all');
subsel = input('?');
switch subsel
    case 1
        disp('Subjects')
        disp(subj_all')
        subsel1 = input('enter subject number?');
end
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
clear datafile2
switch subsel
    case 1
        datafile2{1} = datafile1{subsel1}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    case 2
        datafile2 = datafile1;
end

%% Preprocess
for i = 1:size(datafile2,1)
    
    datafile = datafile2{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile, '/');
    subj = datafile(Index(6)+1:Index(7)-1);
    Index = strfind(datafile, 'run');
    run  = datafile(Index+3);
    disp(datafile)
    disp(['subj:',subj])
    disp(['Run:',run])
    disp('============');
    
    %-elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    disp('============');
    
    outd.sub = fullfile(outdir,'process_ft',subj, tag, run);
    if exist(outd.sub, 'file') == 0
        mkdir(outd.sub);   %create the directory
    end
    cd(outd.sub)
    disp(['outputdir:',outd.sub])
    disp('============');

    %%
    ic_selection = 1; % 1: manual, 2: automated
    Run_preprocess
    
    
    
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













%%

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
flag.speech = 2;

%% Initial settings
% set(0,'DefaultFigureWindowStyle','docked');
% set(0,'DefaultFigureWindowStyle' , 'normal')
% set(gcf,'units','points','position',[500,500,500,500]);

cd '/data/MEG/Research/Aqil_Izadi/Scripts';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/data/MEG/Research/Aqil_Izadi/Raw_data';
%- Output dir
outdir = '/data/MEG/Research/Aqil_Izadi/process';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW_MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
disp('1: SD')
disp('2: Picture naming');
disp('3: Story-math');
task = input('Enter the task: ');
switch task
    case 1
        %- Semantic decision
        tag = 'SD';
    case 2
        %- Visual picture naming
        tag = 'PN';
    case 3
        %- Story math
        tag = 'SM';
end

%% analysis flag
flag.freq = 1;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis

%% Listing data
clc
% if exist(['datalog_',tag,'.mat'], 'file') == 2
%     load(['datalog_',tag,'.mat'])
% else
% d = rdir([datadir,['/**/',ytag,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
% Per year
clear datafolder datafile subj_all
datafile1 = [];
d = rdir([indir,['/**/tSSS/*',tag,'*_raw.fif']]);
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    subj_all{i} = [num2str(i), ': ', datafile{i}(Index(6)+1:Index(7)-1)];
end
datafile1 = vertcat(datafile1,datafile);
datafile1 = datafile1';

%     save(['datalog_',tag,'.mat'],'datafile1','subj_all')
% end
disp(datafile1)
disp('============');

%%
disp('1: choose specific subject');
disp('2: do all');
subsel = input('?');
switch subsel
    case 1
        disp('Subjects')
        disp(subj_all')
        subsel1 = input('enter subject number?');
end
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
clear datafile2
switch subsel
    case 1
        datafile2{1} = datafile1{subsel1}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    case 2
        datafile2 = datafile1;
end

%%
for i = 1:size(datafile2,1)
    
    datafile = datafile2{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile, '/');
    subj = datafile(Index(6)+1:Index(7)-1);
    Index = strfind(datafile, 'run');
    run  = datafile(Index+3);
    disp(datafile)
    disp(['subj:',subj])
    disp(['Run:',run])
    disp('============');
    
    %-elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    disp('============');
    
    %%
    %     if year>=4
    %         yttag = 'older';
    %     else
    %         yttag = ytag{1};
    %     end
    outd.sub = fullfile(outdir,'process_ft',subj, tag, run);
    if exist(outd.sub, 'file') == 0
        mkdir(outd.sub);   %create the directory
    end
    cd(outd.sub)
    disp(['outputdir:',outd.sub])
    disp('============');
    
    %% Preprocesssing
    ic_selection = 1; % 1: manual, 2: automated
    switch task
        case {1,2}
            Run_preprocess
        case 3
            Run_preprocess_motor
            disp('1: left');
            disp('2: right');
            side = input('?');
            switch side
                case 1
                    cln_data = cln.left;
                case 2
                    cln_data = cln.right;
            end
    end
    
    %% Speech
    if flag.speech ==1
        Run_speech
    end
    
    %% Notch filtering of 30Hz
    if flag.notch ==1
        Run_notch
    end
    
    %% FFT & TFR
    fmax = 40;
    datain = cln_data;
    if flag.freq == 1
        stag = 'tsk_baseline';
        Run_freq
        disp(['time_of_interest:',num2str(time_of_interest),'sec']);
        disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
        L = 0.3;
    end
    
    %% Time-locked
    disp(['suggested time interval:',num2str(time_of_interest), '+-', num2str(L),' Sec']);
    toi = [-0.3,0; 0.3,2];
    ep_data = vy_epoch(datain, toi);
    
    cfg = [];
    ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);
    
    %% Grand Mean
    if flag.gave == 1
        Run_grandmean
    end
    
    %% Source analysis
    outputmridir = fullfile(outdir,'process_ft', subj,'anat'); % output dir
    if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
    
    %% Processing
    disp('1: Surface-based')
    disp('2: Volumetric');
    %     analysis = input('Eneter the analysis: ');
    analysis = 2;
    
    switch analysis
        case 1
            mtag = 'source_surf';
            disp('1: Surface-based SPM source analysis');
            disp('2: Export ft to Brainstorm');
            disp('3: Surface-based bf');
            method = input('Method: ');
        case 2
            % end
            disp('1: LCMV')
            disp('2: Network+Connectvity');
            disp('3: DICS Source');
            %             method = input('Method: ');
            method = 3;
            %-
            %             disp('1: Low-res grid')
            %             disp('2: High-res grid')
            %             meshgridres = input('Mesh grid: ');
            meshgridres = 1;
            
            rl = 2;
    end
    
    method = 3;
    switch method
        case 1
            % Surface-based analysis
            Run_surfacebased_ECP
        case 2
            % Volumetric-based analysis
            anatomy_check_flag = 2;
            flag.meshgrid_sel = 1;
            flag.anatomy_check = 1;
            flag.savetag=2;
            %             Run_volumetric
            Run_volumetric_ECP
        case 3
            anatomy_check_flag = 2;
            flag.meshgrid_sel = 1;
            flag.anatomy_check = 1;
            flag.savetag=2;
            % Conn analysis
            Run_anatomy_SD
            Run_sourceanalysis
            Run_extractvirtualsensor
            Run_connanalysis
            Run_networkanalysis %- network analysis (network degrees)
            Run_vis_netconn_analysis %- mapping
            %             Run_dicsanalysis
    end
    disp('============');
    
    %%
    %     pause
    close all
    disp([datafile,' ,was completed'])
    
end

