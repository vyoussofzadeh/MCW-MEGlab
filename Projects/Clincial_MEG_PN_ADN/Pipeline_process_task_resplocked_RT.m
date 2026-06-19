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

%% Initial settings
set(0,'DefaultFigureWindowStyle','docked');
% set(gcf,'units','points','position',[500,500,500,500]);

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

%%
disp('1: Definition naming')
disp('2: Picture naming');
task = input('Eneter the task: ');

switch task
    case 1
        % - Auditory definition naming
        tag = 'DFN'; tag1 = 'dfn';
        Evnt_IDs = 1; % questions
        disp('1: 2019')
        disp('2: 2018');
        disp('3: 2017');
        disp('4: 2016')
        disp('5: 2015');
        disp('6: 2014');
        disp('7: 2013')
        disp('8: 2012');
        disp('9: 2011');
        disp('10: all');
        year = input('Year data were acquired: ');
    case 2
        % - Visual picture naming
        tag = 'PN';
        Evnt_IDs = 3; % 3: images, 2: scrambled images
        disp('1: 2019')
        disp('2: 2018');
        disp('3: 2017');
        year = input('Year data were acquired: ');
end

%%
clear ytag;
switch year
    case 1
        stag = {'19'};ytag = {'19'};
    case 2
        stag = {'18'};ytag = {'18'};
    case 3
        stag = {'17'};ytag = {'17'};
    case 4
        stag = {'16'};ytag = {'16'};
    case {5,6}
        stag = {'up'};ytag = {'14_5'};
        %     case 6
        %         stag = {'up'};ytag = {'14'};
    case 7
        stag = {'13'};ytag = {'13'};
    case 8
        stag = {'12'};ytag = {'12'};
    case 9
        stag = {'11'};ytag = {'11'};
    case 10
        stag = {'all'};ytag = {'all'};
end
disp('============');

%% Listing data
% d = rdir([datadir,['/**/',ytag,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
% Per year
clear datafolder datafile subj_all
datafile1 = [];
for j=1:numel(stag)
    stag1 = stag{1,j};
    d = rdir([indir,['/**/',stag1,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
    %     d = rdir([indir,['/**/',ytag1,'*/','sss','/*',tag,'*/*raw_sss.fif']]);
    if exist('tag1','var')
        d1 = rdir([indir,['/**/',stag1,'*/','sss','/*',tag1,'*/*raw_tsss.fif']]); d=[d;d1];
    end
    for i=1:length(d)
        [pathstr, name] = fileparts(d(i).name);
        datafolder{i} = pathstr;
        datafile{i} = d(i).name;
        Index = strfind(datafile{i}, '/');
        subj_all_num{i} = [num2str(i), ': ', datafile{i}(Index(5)+1:Index(6)-1)];
        subj_all{i} = datafile{i}(Index(5)+1:Index(6)-1);
    end
    datafile1 = vertcat(datafile1,datafile);
end
datafile1 = datafile1';
disp(datafile1)
disp('============');

%% Finding dublicates (runs 1-2)
[subj_all_uniq, ~, uni_idx] = unique(subj_all,'first');
for i=1:length(subj_all_uniq)
    subj_all_uniq_num{i} = [num2str(i), ': ', subj_all_uniq{i}];
end

%%
disp('1: choose specific subject(s)');
disp('2: do all');
subsel = input('?');
switch subsel
    case 1
        disp('Subjects')
        %         disp(datafile1)
        disp(subj_all_uniq_num')
        subsel1 = input('enter subject number?');
        nsub = length(subsel);
    case 2
        nsub = length(subj_all_uniq);
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

% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 2])
% print('test','-depsc', '-r100')

%%
clear datafile2
switch subsel
    case 1
        if subsel1 > 1
            idx = find(uni_idx==subsel1); nrun = length(idx);
            for i=1:nrun
                datafile2{i,:} = datafile1{idx(i)}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
            end
        else
            datafile2{1} = datafile1{subsel1}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
            nrun = 1;
        end
        disp(datafile2)
    case 2
        datafile2 = datafile1;
        disp(subj_all'); nrun = 1;
end

%%
for i = 1:nsub
    
    close all
    datafile = datafile2{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile, '/');
    subj = datafile(Index(5)+1:Index(6)-1);
    Date  = datafile(Index(6)+1:Index(7)-1);
    disp('============');
    disp([num2str(i),'/',num2str(nsub)])
    disp(datafile)
    disp(['subj:',subj])
    disp(['Date:',Date])
    disp('============');
    
    %-elec/grad
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    disp('============');
    
    %%
    %     if year>=4
    %         yttag = 'older';
    %     else
    yttag = ytag{1};
    %     end
    outd.sub = fullfile(outdir,'ft_process',yttag, subj, tag);
    if exist(outd.sub, 'file') == 0
        mkdir(outd.sub);   %create the directory
    end
    cd(outd.sub)
    %     pause(2)
    disp(['outputdir:',outd.sub])
    disp('============');
    
    %% Preprocesssing
    clear cln_data
    if nrun > 1
        datafile = datafile2{end};
    end
    ic_selection = 1; % 1: manual, 2: automated
    Run_preprocess
    
    %% Speech
    Run_speech
    
    disp([num2str(mean(tt)),'+-', num2str(std(tt)),' Sec'])
%     length(ep_data.pst.sampleinfo)/length(f_data.sampleinfo)
    
    %% Trial-by-trail variability correction based on speech onset of good trials (response-lock analysis)
    %     Run_TBT
    
end

