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
addpath('./data_full')
addpath('./run')

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/helper')

%% Read raw data
cd(indir)
tag = 'SD'; %- Semantic decision

clc
if exist(['datalog_ic_',tag,'.mat'], 'file') == 2
    load(['datalog_ic_',tag,'.mat'])
else
    clear datafolder datafile subj_all sub
    datafile_fif = [];
    d = rdir([indir,['/**/ICA/*',tag,'*/ec*_ica_clean.fif']]);
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
    save(['datalog_ic_',tag,'.mat'],'datafile_fif','subj_all', 'sub')
end
disp(datafile_fif)

%% BS
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3';

addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

BS_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/';
BS_data_dir = fullfile(BS_dir,'data_full');
protocol = fullfile(BS_dir, 'data_full/protocol.mat');

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
    datadir_sub = fullfile(BS_dir,'data_full/',['EC',sub_sel]);
    cd(datadir_sub)
    idx = strfind(datadir_sub,'/');
    okrun = find(contains(no_anat,datadir_sub(idx(end)+1:end))==1);
    if  ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw']) ...
            && isempty(okrun) && ...
            ~isfolder(['@rawEC',sub_sel, '_SD_', run_sel, '_raw'])...
            && ~isfolder(['@rawEC',sub_sel, '_SD_', run_sel, '_raw_ica_clean']) ...
            && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw_ica_clean'])...
            && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw_tsss_ica_clean'])...
            && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw_tsss_ica_clean_low_clean'])...
            && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw_tsss']) ...
            && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_elecfix_raw_ica_clean'])
        iSubject = find(contains(unq_bs_subj, sub_sel)==1);
        RawFiles = datafile_fif{ii};
        disp(RawFiles)
        pause,
%         OutputFiles = import_raw(RawFiles, 'FIF', iSubject, [], []); % GUI-based
        ecp_import_raw(RawFiles, iSubject); % no GUI
        
    end
end

%% Preprocess raw data
% Unload('fieldtrip');

d1 = rdir(fullfile(BS_data_dir,'/EC*/@raw*/*_raw_ica_clean.mat'));
d2 = rdir(fullfile(BS_data_dir,'/EC*/@raw*/*_raw_tsss_ica_clean.mat'));
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
    idx2 = strfind(name,'run');
    if isempty(idx2)
        idx2 = strfind(name,'Run');
    end
    sub_sel = name(idx+2:idx+5);
    run_sel = name(idx2+3);
    cd(pathstr2)
    if  ~isfolder(['@rawec',sub_sel, '_SD_run', run_sel, '_raw_ica_clean_low_clean']) ...
            && ~isfolder(['@rawec',sub_sel, '_SD_Run', run_sel, '_raw_ica_clean_low_clean']) ...
            && ~isfolder(['@rawec',sub_sel, '_SD_run', run_sel, '_raw_ica_clean_low']) ...
            && ~isfolder(['@rawec',sub_sel, '_SD_run', run_sel, '_raw_tsss_ica_clean_low_clean']) ...
            && ~isfolder(['@rawec',sub_sel, '_SD_run', run_sel, '_elecfix_raw_ica_clean_low'])
%         ...
%             && ~isfolder(['@rawec',sub_sel, '_SD_Run', run_sel, '_raw_ica_clean'])
        cd(pathstr2)
        disp(sFiles)
        disp(i)
        pause,
        bs_preprocess(sFiles)
    end
end

%% Import epoched data
% Load('fieldtrip');

% db_reload_database('current',1)
d1 = rdir(fullfile(BS_data_dir,'/**/@raw*/*clean_low_clean.mat'));
d2 = rdir(fullfile(BS_data_dir,'/**/@raw*/*raw_ica_clean_low.mat'));

% d2 = rdir(fullfile(BS_data_dir,'/**/@raw*/*tsss_ica_clean.mat'));

d = [d1; d2];
% d = d1;
% d = rdir(fullfile(BS_data_dir,'/**/@raw*/*raw_ica_clean_low.mat'));

% Exception: ec1134 (no eog?)

clc
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    [pathstr2, name2] = fileparts(pathstr);
    [pathstr3, name3] = fileparts(pathstr2);
    sFiles = {fullfile(name3, name2, [name, '.mat'])};
    idx = strfind(name,'ec'); idx1 = strfind(name,'_');
    id2 = strfind(name,'_run');
    if isempty(id2)
        id2 = strfind(name,'_Run');
    end
    sub_sel = name(idx+2:idx+5);
    run_sel = name(id2+1:id2+4);
    iSubject = find(contains(unq_bs_subj, sub_sel)==1);
    cd(pathstr2)
    if  ~isfolder(['ec',sub_sel, '_SD_', run_sel, '_raw_ica_clean_low_clean']) ...
            && ~isfolder(['ec',sub_sel, '_SD_', run_sel, '_raw_tsss_ica_clean_low_clean'])
        disp(sFiles)
        
        pause,
        %         import_raw_to_db(sFiles{1}); % GUI-based
        
        % Define epoching parameters
        EpochTime = [-0.5, 2];  % Epoch from -100 ms to +300 ms around the event
        subjectName = ['EC',sub_sel];
        
        % Create epochs
        sFiles_3 = bst_process('CallProcess', 'process_import_data_event', sFiles{1}, [], ...
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
        
        % Reject bad trials
        bst_process('CallProcess', 'process_ft_reject_trials_edit', sFiles_3, [], 'sensortype', 'MEG', 'mrej', 'Auto', 'arej', 0.9); % Process: FieldTrip: process_ft_reject_trials kurtosis > 15
        
        
        % Create epochs
        sFiles_2 = bst_process('CallProcess', 'process_import_data_event', sFiles{1}, [], ...
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
        
        % Reject bad trials
        bst_process('CallProcess', 'process_ft_reject_trials_edit', sFiles_2, [], 'sensortype', 'MEG', 'mrej', 'Auto', 'arej', 0.9); % Process: FieldTrip: process_ft_reject_trials kurtosis > 15
        
    end
end

%% Refine the headpoint
% 1) Set your root BS "data_full" directory
rootDir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_full';

% 2) List all subfolders named "ECxxxx"
dSubs = dir(fullfile(rootDir, 'EC*')); % e.g. "EC1002", "EC1077", etc.
dSubs = dSubs([dSubs.isdir]);         % keep only directories

% 3) For each subject subfolder, gather channel_vectorview306_acc1.mat files
channelFiles = {};   % we'll store them in a cell array
for i = 1:length(dSubs)
    subjFolder = fullfile(rootDir, dSubs(i).name); % e.g. /data_full/EC1002
    
    % Find all channel_vectorview306_acc1.mat files inside this subject folder
    % (Sometimes they may be in sub-subfolders like @raw etc.)
    pattern   = fullfile(subjFolder, '/ec*/', 'channel_vectorview306_acc1.mat');
    foundList = dir(pattern);
    for j = 1:length(foundList)
        channelFiles{end+1} = fullfile(foundList(j).folder, foundList(j).name);
    end
end

fprintf('Found %d channel_vectorview306_acc1.mat files total.\n', numel(channelFiles));

% 4) Refine each channel file with channel_align_auto()
ChannelMat = [];
isWarning  = 1; % show warnings in console
isConfirm  = 0; % no user prompts
tolerance  = 0;
ignorepercent = '0';

for iFile = 1:numel(channelFiles)
    ChannelFile = channelFiles{iFile};
    fprintf('\n[%d/%d] Refining alignment for:\n  %s\n', iFile, numel(channelFiles), ChannelFile);
    
    % Call Brainstorm function (make sure Brainstorm is on path and BST is open/started)
    channel_align_auto_ecp(ChannelFile, ChannelMat, isWarning, isConfirm, tolerance, ignorepercent);
    
    close all
    
    % You can optionally reload BST database after each run:
    % db_reload_database('current',1);
end

fprintf('\nDone! Processed %d channel files.\n', numel(channelFiles));


%% Est. head model
overWriteHM = true;

d1 = rdir(fullfile(BS_data_dir,'/*/ec*_ica_clean_low_clean/channel_vectorview306_acc1.mat'));
d2 = rdir(fullfile(BS_data_dir,'/*/ec*_raw_ica_clean_low_clean/channel_vectorview306_acc1.mat'));
d3 = rdir(fullfile(BS_data_dir,'/*/EC*raw_ica_low_clean/channel_vectorview306_acc1.mat'));
d4 = rdir(fullfile(BS_data_dir,'/*/ec*raw_ica_clean_low_clean/channel_vectorview306_acc1.mat'));
d5 = rdir(fullfile(BS_data_dir,'/*/ec*elecfix_raw_ica_clean*/channel_vectorview306_acc1.mat'));
d6 = rdir(fullfile(BS_data_dir,'/*/ec*_raw_ica_clean_low/channel_vectorview306_acc1.mat'));
d7 = rdir(fullfile(BS_data_dir,'/*/ec*_raw_ica_clean_low_clean/channel_vectorview306_acc1.mat'));
d8 = rdir(fullfile(BS_data_dir,'/*/ec*_raw_ica_clean_low/channel_vectorview306_acc1.mat'));

d = [d1;d2;d3;d4;d5;d7;d8];

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
    headmodelFile = 'headmodel_surf_os_meg.mat';
    cd(a)
    if ~~exist(fullfile(BS_dir,'anat',OPTIONS.CortexFile),'file') && ~exist(headmodelFile,'file') || overWriteHM == true
        if overWriteHM == true
            delete(headmodelFile);
            pause(2)
        end
%         pause
        bst_headmodeler(OPTIONS);
    end
    disp(a)
end
% % db_reload_database('current',1)



