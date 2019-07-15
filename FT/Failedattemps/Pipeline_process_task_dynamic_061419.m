clear; clc, close('all'); warning off

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
        tag = 'DFN'; tag1 = 'dfn';
        Evnt_IDs = 1; % questions
    case 2
        % - Visual picture naming
        tag = 'PN';
        Evnt_IDs = 3; % 3: images, 2: scrambled images
end

%% analysis flag
flag.freq = 0;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis

%%
disp('1: 2019')
disp('2: 2018');
disp('3: 2017');
disp('4: 2016')
disp('5: 2015');
disp('6: 2014');
disp('7: 2013')
disp('8: 2012');
disp('9: 2011');
disp('10: older');
year = input('Year data were acquired: ');

clear ytag;
switch year
    case 1
        ytag = {'19'};
    case 2
        ytag = {'18'};
    case 3
        ytag = {'17'};
    case 4
        ytag = {'16'};
    case 5
        ytag = {'up'};
    case 6
        ytag = {'up'};
    case 7
        ytag = {'13'};
    case 8
        ytag = {'12'};
    case 9
        ytag = {'11'};
end
disp('============');

%% Listing data
% d = rdir([datadir,['/**/',ytag,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
% Per year
clear datafolder datafile subj_all
datafile1 = [];
for j=1:numel(ytag)
    ytag1 = ytag{1,j};
    d = rdir([indir,['/**/',ytag1,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
    %     d = rdir([indir,['/**/',ytag1,'*/','sss','/*',tag,'*/*raw_sss.fif']]);
    if exist('tag1','var')
        d1 = rdir([indir,['/**/',ytag1,'*/','sss','/*',tag1,'*/*raw_tsss.fif']]); d=[d;d1];
    end
    for i=1:length(d)
        [pathstr, name] = fileparts(d(i).name);
        datafolder{i} = pathstr;
        datafile{i} = d(i).name;
        Index = strfind(datafile{i}, '/');
        subj_all{i} = [num2str(i), ': ', datafile{i}(Index(5)+1:Index(6)-1)];
    end
    datafile1 = vertcat(datafile1,datafile);
end
datafile1 = datafile1';
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
    subj = datafile(Index(5)+1:Index(6)-1);
    Date  = datafile(Index(6)+1:Index(7)-1);
    disp(datafile)
    disp(['subj:',subj])
    disp(['Date:',Date])
    disp('============');
    
    %-elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    disp('============');
    
    %%
    if year>=4
        yttag = 'older';
    else
        yttag = ytag{1};
    end
    outd.sub = fullfile(outdir,'ft_process',yttag, subj, tag);
    if exist(outd.sub, 'file') == 0
        mkdir(outd.sub);   %create the directory
    end
    cd(outd.sub)
    disp(['outputdir:',outd.sub])
    disp('============');
    
    %% Preprocesssing
    Run_preprocess
    %%
    cln_data1 = cln_data;
    cln_data = cln_data1;
    
    cfg = [];
    cfg.savefile = [];
    cfg.saveflag = 2;
    cfg.foilim = [2 40];
    cfg.plotflag  = 1;
    vy_fft(cfg, cln_data);
    grid on
    grid minor
    title('Before band-stop filtering');
    
    %% 
    fsb = input('Enter the sop-band frequency?');
    cfg = [];
    cfg.bsfilter = 'yes';
%     cfg.bsfreq = [29 32]; % or whatever you deem appropriate
    cfg.bsfreq = [fsb-1 fsb+1]; % or whatever you deem appropriate
%     cfg.bsfreq = [8 12;29 32]; % or whatever you deem appropriate
    cln_data = ft_preprocessing(cfg, cln_data);
%     cfg.bsfreq = [2 12]; % or whatever you deem appropriate


    cfg = [];
    cfg.savefile = [];
    cfg.saveflag = 2;
    cfg.foilim = [2 40];
    cfg.plotflag  = 1;
    vy_fft(cfg, cln_data);
    grid on
    grid minor
    title('After band-stop filtering')
    
    %% FFT & TFR
    if flag.freq == 1
        Run_freq
    end
    
    %% Timelock
    if flag.time == 1
        switch task
            case 1
                %                 toi = [-0.3,0;1.1,1.8]; % Best of DN
                toi = [-0.3,0;1.2,1.6]; % Best of DN
                %                 toi = [-0.3,0;0.8,1.2]; % Best of DN
                %                 toi = [-0.3,0;0.4,0.8]; % Best of DN
                %                 toi = [-0.3,0;1.5,2]; % Best of DN
            case 2
                %                 toi = [0,0.3;0.6,1.2]; % Best for PN, left IFG
%                 toi = [0,0.3;0.4,1.2]; % Best for PN, left IFG
%                 toi = [-0.3,0;0.4,1.2]; % Best for PN, left IFG
                toi = [-0.3,0;0.4,1.2]; % Best for PN, left IFG
%                 toi = [-0.3,0;0.4,0.8]; % Best for PN, left IFG
                %                 toi = [0,0.3;0.4,1.6];
%                 toi = [0,0.3;0.6,1.2];
        end
        Run_time
    end
    
    %% Grand Mean
%     if flag.gave == 1
%         Run_grandmean
%     end
    %% Source analysis
    outputmridir = fullfile(outdir,'ft_process',yttag, subj,'anat'); % output dir
    if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
    
    %%
    clear toi   
    switch task
        case 1
            toi{1} = [-0.4,0;0.4,0.8];
            toi{2} = [-0.4,0;0.8,1.2];
%             toi{3} = [-0.4,0;1.2,1.6];
            %     toi{4} = [-0.4,0;1.6,2];
        case 2
%             toi{1} = [0,0.3;0.4,0.8];
            toi{1} = [0,0.3;0.6,1.2];
%             toi{3} = [0,0.3;1.2,1.6];
            %     toi{4} = [0,0.4;1.6,2];
    end
    
    %% Processing
    disp('1: Surface-based')
    disp('2: Volumetric');
    analysis = input('Eneter the analysis: ');
    
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
            %-
%             disp('1: Low-res grid')
%             disp('2: High-res grid')
%             meshgrid = input('Mesh grid: ');
            meshgrid = 1;
    end
    
    switch analysis
        case 1
            % Surface-based analysis
            Run_surfacebased
            
        case 2
            
            % Volumetric-based analysis
            anatomy_check_flag = 2;
            Run_volumetric
    end
    disp('============');
    
    %%
    pause
    close all
    disp([datafile,' ,was completed'])
    
end




