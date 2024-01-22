%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (preprocessing)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 01/19/2024

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
cd('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD')
addpath('./run')
Run_setpath
addpath('./data')
addpath('./run')

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/helper')

%% Loading up raw data
cd(indir)
tag = 'SD'; %- Semantic decision

clc
if exist(['datalog_',tag,'.mat'], 'file') == 2
    load(['datalog_',tag,'.mat'])
else
    clear datafolder datafile subj_all sub
    datafile_fif = [];
    d1 = rdir([indir,['/**/tSSS/*',tag,'*_raw.fif']]);
    d2 = rdir([indir,['/**/tSSS/*',tag,'*_raw_tsss.fif']]);
    d = [d1;d2];
    for i=1:length(d)
        [pathstr, name] = fileparts(d(i).name);
        datafolder{i} = pathstr;
        datafile{i} = d(i).name;
        Index = strfind(datafile{i}, '/');
        sub{i} = datafile{i}(Index(6)+1:Index(7)-1);
        subj_all{i} = [num2str(i), ': ', datafile{i}(Index(6)+1:Index(7)-1)];
    end
    datafile_fif = vertcat(datafile_fif,datafile);
    datafile_fif = datafile_fif';
    save(['datalog_',tag,'.mat'],'datafile_fif','subj_all', 'sub')
end
disp(datafile_fif)

%% BS
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3';

addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

BS_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/';
BS_data_dir = fullfile(BS_dir,'data_all_subjects');
protocol = fullfile(BS_dir, 'data/protocol.mat');

%%
db_reload_database('current',1)
load(protocol);
Subj_bs = ProtocolSubjects.Subject;

L = length(Subj_bs);
k = 1;
clear subjs_bs
for i=1:length(Subj_bs)
    if ~contains(Subj_bs(i).Name, 'Group_analysis')
        datafile{k} = Subj_bs(i).FileName;
        subjs_bs{k} = Subj_bs(i).Name;
        k=1+k;
    end
end
unq_bs_subj = unique(subjs_bs);

%%
no_anat = {'EC1036'
    'EC1037'
    'EC1038'
    'EC1040'
    'EC1045'
    'EC1049'
    'EC1061'
    'EC1065'
    'EC1085'
    'EC1094'
    'EC1096'
    'EC1110'
    'EC1111'
    'EC1090'
    };

temp_anat = {'EC1090'};

sub_all1 = setdiff(unq_bs_subj,no_anat);

%% Import raw data
clc
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/BS_func_edits')

L_data =length(datafile_fif);
not_imported = [];
for ii=1:L_data
    [~, name] = fileparts(datafile_fif{ii});
    idx = strfind(name,'_'); sub_sel = name(3:idx(1)-1); run_sel = name(idx(2)+1:idx(3)-1);
    datadir_sub = fullfile(BS_dir,'data/',['EC',sub_sel]);
    cd(datadir_sub)
    idx = strfind(datadir_sub,'/');
    okrun = find(contains(no_anat,datadir_sub(idx(end)+1:end))==1);
    if  ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw']) ...
            && isempty(okrun) && ...
            ~isfolder(['@rawEC',sub_sel, '_SD_', run_sel, '_raw'])...
            && ~isfolder(['@rawEC',sub_sel, '_SD_', run_sel, '_raw_ica_clean']) ...
            && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw_ica_clean'])...
            && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw_tsss']) ...
            && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_elecfix_raw_ica_clean'])
        iSubject = find(contains(unq_bs_subj, sub_sel)==1);
        RawFiles = datafile_fif{ii};
        disp(RawFiles)
%         pause,
%         OutputFiles = import_raw(RawFiles, 'FIF', iSubject, [], []); % GUI-based
        ecp_import_raw(RawFiles, iSubject); % no GUI
        
    end
end

%% Preprocess raw data
d1 = rdir(fullfile(BS_data_dir,'/EC*/@raw*/*_raw.mat'));
d2 = rdir(fullfile(BS_data_dir,'/EC*/@raw*/*_raw_tsss.mat'));
d = [d1;d2];
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    [pathstr2, name2] = fileparts(pathstr);
    [pathstr3, name3] = fileparts(pathstr2);
    sFiles = {fullfile(name3, name2, [name, '.mat'])};
    idx = strfind(name,'ec');
    if isempty(idx)
        idx = strfind(name,'EC');
    end
    %     idx1 = strfind(name,'_');
    idx2 = strfind(name,'run');
    if isempty(idx2)
        idx2 = strfind(name,'Run');
    end
    sub_sel = name(idx+2:idx+5);
    run_sel = name(idx2+3);
    cd(pathstr2)
    if  ~isfolder(['@rawec',sub_sel, '_SD_run', run_sel, '_raw_low']) ...
            && ~isfolder(['@rawEC',sub_sel, '_SD_run', run_sel, '_raw_low'])...
            && ~isfolder(['@rawEC',sub_sel, '_SD_run', run_sel, '_elecfix_raw_low'])...
            && ~isfolder(['@rawec',sub_sel, '_SD_run', run_sel, '_elecfix_raw_low'])...
            && ~isfolder(['@rawec',sub_sel, '_SD_Run', run_sel, '_raw_low'])...
            && ~isfolder(['@rawec',sub_sel, '_SD_Run', run_sel, '_raw_tsss_low'])...
            && ~isfolder(['@rawec',sub_sel, '_SD_run', run_sel, '_raw_tsss_low'])
        disp(sFiles)
        %         pause,
        %         cd(BS_data_dir)
        bs_preprocess(sFiles)
    end
end

%% Import epoched data
db_reload_database('current',1)
d = rdir(fullfile(BS_data_dir,'/**/@raw*/*_low_clean.mat'));
% d = rdir(fullfile(BS_data_dir,'/**/@raw*/*_low_02_clean.mat'));

clc
for i=1:1%length(d)
    [pathstr, name] = fileparts(d(i).name);
    [pathstr2, name2] = fileparts(pathstr);
    [pathstr3, name3] = fileparts(pathstr2);
    sFiles = {fullfile(name3, name2, [name, '.mat'])};
    idx = strfind(name,'ec'); idx1 = strfind(name,'_');
    sub_sel = name(idx+2:idx+5);
    run_sel = name(idx1(end-5)+1:idx1(end-4)-1);
    iSubject = find(contains(unq_bs_subj, sub_sel)==1);
    cd(pathstr2)
    if  ~isfolder(['ec',sub_sel, '_SD_', run_sel, '_raw_ica_clean_low_clean'])
        disp(sFiles)
        
        %         pause,
        %         import_raw_to_db(sFiles{1}); % GUI-based
        
        % Define epoching parameters
        EpochTime = [-0.3, 2];  % Epoch from -100 ms to +300 ms around the event
        subjectName = ['EC',sub_sel];
        
        % Create epochs
        bst_process('CallProcess', 'process_import_data_event', sFiles{1}, [], ...
            'subjectname', subjectName, ...
            'condition', '', ...
            'eventname', '3', ...
            'timewindow', [], ...
            'epochtime', EpochTime, ...
            'createcond', 0, ...
            'ignoreshort', 1, ...
            'usectfcomp', 1, ...
            'usessp', 1, ...
            'freq', 1000, ...
            'baseline', []);
        
        % Create epochs
        bst_process('CallProcess', 'process_import_data_event', sFiles{1}, [], ...
            'subjectname', subjectName, ...
            'condition', '', ...
            'eventname', '2', ...
            'timewindow', [], ...
            'epochtime', EpochTime, ...
            'createcond', 0, ...
            'ignoreshort', 1, ...
            'usectfcomp', 1, ...
            'usessp', 1, ...
            'freq', 1000, ...
            'baseline', []);
    end
end

%% Est. head model
d1 = rdir(fullfile(BS_data_dir,'/*/ec*_ica_clean_low_clean/channel_vectorview306_acc1.mat'));
d2 = rdir(fullfile(BS_data_dir,'/*/ec*_raw_ica_clean_low_clean/channel_vectorview306_acc1.mat'));
% d2 = rdir(fullfile(BS_data_dir,'/*/ec*_low_02_clean/channel_vectorview306_acc1.mat'));
d3 = rdir(fullfile(BS_data_dir,'/*/EC*raw_ica_low_clean/channel_vectorview306_acc1.mat'));
d4 = rdir(fullfile(BS_data_dir,'/*/ec*raw_ica_clean_low_clean/channel_vectorview306_acc1.mat'));
d5 = rdir(fullfile(BS_data_dir,'/*/ec*elecfix_raw_ica_clean*/channel_vectorview306_acc1.mat'));

d = [d1;d2;d3; d4; d5];

OPTIONS = [];
OPTIONS.comment = 'Overlapping spheres';
OPTIONS.MEGMethod =  'os_meg';
OPTIONS.EEGMethod  ='';
OPTIONS.ECOGMethod = '';
OPTIONS.SEEGMethod = '';
OPTIONS.SaveFile = 1;

for ii=1:length(d)
    [a, ~] = fileparts(d(ii).name);
    [c,~] = fileparts(a);
    [e,f] = fileparts(c);
    OPTIONS.HeadModelFile =  a;
    OPTIONS.HeadModelType  = 'surface';
    A = load(d(ii).name);
    OPTIONS.Channel = A.Channel;
    OPTIONS.CortexFile = fullfile(f,'tess_cortex_pial_low.mat');
    OPTIONS.HeadFile =  fullfile(f,'tess_head_mask.mat');
    cd(a)
    if ~~exist(fullfile(BS_dir,'anat',OPTIONS.CortexFile),'file') && ~exist('headmodel_surf_os_meg.mat','file')
        bst_headmodeler(OPTIONS);
    end
    disp(a)
end
db_reload_database('current',1)
