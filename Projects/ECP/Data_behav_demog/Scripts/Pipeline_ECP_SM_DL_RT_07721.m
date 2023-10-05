clear; clc, close('all'); warning off

%% Flags
flag.preprocessing.filtering = 1;
flag.trialinfo.filtering = 1;
flag.preprocessing.artifact = 1;
flag.preprocessing.ica = 1;
flag.notch = 0;
flag.freq = 0;     % TFR & FFT
flag.toi = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis
flag.speech = 0;             % speech analysis
flag.sourceanalysis = 1;     % source analysis
flag.anatomy = 1;            % grand average analysis
flag.anatomy_check = 0;
flag.meshgrid_sel = 1; % high-res = 2, low-res = 1;
flag.overwrite = 1;

%% Initial settings
% set(0,'DefaultFigureWindowStyle','normal');
% set(gcf,'units','points','position',[500,500,500,500]);
% set(0, 'DefaultFigureRenderer', 'painters');
% set(gcf, 'renderer', 'painters');
% set(0, 'DefaultFigureRenderer', 'opengl');

cd '/data/MEG/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));
rmpath('./Failedattemps');

indir = '/group/jbinder/ECP/MEG/MEG_Work';
%- Output dir
outdir = '/data/MEG/Research/ECP/';

ECP_scriptdir = '/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP';
% ECP_datadir = '/data/MEG/Research/ECP/';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
%- Picture naming
tag = [];
tag.task = 'SM';

%% Listing data
clc
if exist(['datalog_',tag.task,'.mat'], 'file') == 2
    load(['datalog_',tag.task,'.mat'])
else
% d = rdir([datadir,['/**/',ytag,'*/','sss','/*',tag,'*/*raw_tsss.fif']]);
% Per year
clear datafolder datafile subj_all
datafile1 = [];
d = rdir([indir,['/**/tSSS/*',tag.task,'*_raw.fif']]);
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    subj_all{i} = [num2str(i), ': ', datafile{i}(Index(6)+1:Index(7)-1)];
end
datafile1 = vertcat(datafile1,datafile);
datafile1 = datafile1';

save(['datalog_',tag.task,'.mat'],'datafile1','subj_all')
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

%%
for dd = 1:size(datafile2,1)
    
    datafile = datafile2{dd}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
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
    outd.sub = fullfile(outdir,'MEG_work_ft',subj, tag.task, run);
    if exist(outd.sub, 'file') == 0
        mkdir(outd.sub);   %create the directory
    end
    cd(outd.sub)
    disp(['outputdir:',outd.sub])
    disp('============');
    
    %%
    %     %% Preprocesssing
    %     ic_selection = 1; % 1: manual, 2: automated
    %     Run_preprocess_SM_3
    Run_tringger_SM
    
    %%
    idx_str = trlInfoColDescr(:,2)==1;
    idx_math = trlInfoColDescr(:,2)==2;
    
    %%
    %     str_run = trlInfoColDescr(idx_str,3); str_run_uq = unique(str_run);
    %     for iq = 1:length(str_run_uq)
    %         idx = find(str_run == str_run_uq(iq));
    %         idx_str_run_uq(iq) = idx_str(idx(1));
    %     end
    %     length(idx_str_run_uq)
    %
    %     math_run = trlInfoColDescr(idx_math,3);
    %     math_run_uq = unique(math_run);
    %     for iq = 1:length(math_run_uq)
    %         idx = find(math_run == math_run_uq(iq));
    %         idx_math_run_uq(iq) = idx_math(idx(1));
    %     end
    %     length(idx_math_run_uq)
    
    %% Response time
    savefile = ['RT_',subj,'.mat'];
    if (exist(savefile, 'file') == 2) && (flag.overwrite ~= 1)
        load(savefile);
    else
        
        RT = (trlInfoColDescr(:,17) - trlInfoColDescr(:,15))./2e3;
        RT_str = RT(idx_str); idx_str2 = RT_str > 0; RT_str_uniq = unique(RT_str(idx_str2));
        RT_math = RT(idx_math); idx_math2 = RT_math > 0; RT_math_uniq = unique(RT_math(idx_math2));
        
        m_RT_str_uniq = mean(RT_str_uniq);
        m_RT_math_uniq = mean(RT_math_uniq);
        
        disp(['RT of str was', num2str(m_RT_str_uniq), ' sec']);
        disp(['RT of math was', num2str(m_RT_math_uniq), ' sec']);
        
        save(savefile, 'RT', 'm_RT_str_uniq', 'm_RT_math_uniq');
    end
    
    %% Difficulty level
    savefile = ['DL_',subj,'.mat'];
    if (exist(savefile, 'file') == 2) && (flag.overwrite ~= 1)
        load(savefile);
    else
        %         Run_triggers_SM
        load(['trial_info_',subj,'.mat']);
        idx_str = trlInfoColDescr(:,2)==1;
        idx_math = trlInfoColDescr(:,2)==2;
        
        DL = (trlInfoColDescr(:,7));
        DL_str = DL(idx_str); idx_str2 = find(DL_str > 0); DL_str_uniq = unique(DL_str(idx_str2));
        DL_math = DL(idx_math); idx_math2 = find(DL_math > 0); DL_math_uniq = unique(DL_math(idx_math2));
        
        %         DL_str_uniq
        %         DL_math_uniq
        
        save(savefile, 'DL', 'DL_str_uniq', 'DL_math_uniq');
    end
    
    %% Accuracy
    %     pause
    savefile = ['Acc_',subj,'.mat'];
    if (exist(savefile, 'file') == 2) && (flag.overwrite ~= 1)
        load(savefile);
    else
        load(['trial_info_',subj,'.mat']);
        idx_str = find(trlInfoColDescr(:,2)==1);
        idx_math = find(trlInfoColDescr(:,2)==2);
        
        %         accuracy = length(find(trlInfoColDescr(:,18)==1))/length(trlInfoColDescr)*100;
        accuracy_str  = length(find(trlInfoColDescr(idx_str,18)==1))/length(idx_str)*100;
        accuracy_math = length(find(trlInfoColDescr(idx_math,18)==1))/length(idx_math)*100;
        accuracy = (accuracy_math + accuracy_str)/2;
        
        
        disp(['Accuracy was: %', num2str(accuracy)]);
        disp(['Accuracy of str was: %', num2str(accuracy_str)]);
        disp(['Accuracy of math was: %', num2str(accuracy_math)]);
        
        idx_str_run_uq = [];
        str_run = trlInfoColDescr(idx_str,3); str_run_uq = unique(str_run);
        for iq = 1:length(str_run_uq)
            idx = find(str_run == str_run_uq(iq));
            idx_str_run_uq(iq) = idx_str(idx(1));
        end
        %         length(idx_str_run_uq)
        
        idx_math_run_uq = [];
        math_run = trlInfoColDescr(idx_math,3);
        math_run_uq = unique(math_run);
        for iq = 1:length(math_run_uq)
            idx = find(math_run == math_run_uq(iq));
            idx_math_run_uq(iq) = idx_math(idx(1));
        end
        %         length(idx_math_run_uq)
        
%         idx_str_res = find(trlInfoColDescr(:,21)==15);
%         idx_math_res = find(trlInfoColDescr(:,21)==25);
%         length(find(trlInfoColDescr(idx_str_res,18)==1))/length(idx_str_res)*100
%         length(find(trlInfoColDescr(idx_math_res,18)==1))/length(idx_math_res)*100

        
        accuracy_str  = length(find(trlInfoColDescr(idx_str_run_uq,18)==1))/length(idx_str_run_uq)*100;
        accuracy_math = length(find(trlInfoColDescr(idx_math_run_uq,18)==1))/length(idx_math_run_uq)*100;
        accuracy = (accuracy_math + accuracy_str)/2;
        
        disp(['Accuracy was: %', num2str(accuracy)]);
        disp(['Accuracy of str was: %', num2str(accuracy_str)]);
        disp(['Accuracy of math was: %', num2str(accuracy_math)]);
        %         pause,
        save(savefile, 'accuracy', 'accuracy_str', 'accuracy_math', 'idx_str', 'idx_math');
    end
    close all
end

