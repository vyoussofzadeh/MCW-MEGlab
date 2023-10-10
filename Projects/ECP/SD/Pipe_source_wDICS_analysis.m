%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (preprocessing, source analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 11/20/2022

clear; clc, close('all'); warning off,

%%
analysis_flag.ica = 0;

%% Paths
restoredefaultpath
cd('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD')
addpath('./run')
Run_setpath
addpath('./data')
addpath('./run')
% addpath(genpath('./functions'))

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

%%
% adding BS path
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3';
% bs_path = '/opt/matlab_toolboxes/brainstorm3';

addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

BS_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/';
BS_data_dir = fullfile(BS_dir,'data_all_subjects');
protocol = fullfile(BS_dir, 'data_all_subjects/protocol.mat');

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
L_data =length(datafile_fif);
not_imported = [];
for ii=1:L_data
    [~, name] = fileparts(datafile_fif{ii});
    idx = strfind(name,'_'); sub_sel = name(3:idx(1)-1); run_sel = name(idx(2)+1:idx(3)-1);
    datadir_sub = fullfile(BS_dir,'data_all_subjects/',['EC',sub_sel]);
    cd(datadir_sub)
    idx = strfind(datadir_sub,'/');
    okrun = find(contains(no_anat,datadir_sub(idx(end)+1:end))==1);
    if  ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw']) ...
            && isempty(okrun) && ...
            ~isfolder(['@rawEC',sub_sel, '_SD_', run_sel, '_raw'])...
            && ~isfolder(['@rawEC',sub_sel, '_SD_', run_sel, '_elecfix_raw']) ...
            && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_elecfix_raw'])...
            && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw_tsss'])
        iSubject = find(contains(unq_bs_subj, sub_sel)==1);
        RawFiles = datafile_fif{ii};
        pause(3),
        OutputFiles = import_raw(RawFiles, 'FIF', iSubject, [], []);
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
d1 = rdir(fullfile(BS_data_dir,'/**/@raw*/*raw_low_clean.mat'));
d2 = rdir(fullfile(BS_data_dir,'/**/@raw*/*raw_tsss_low_clean.mat'));
d = [d1;d2];
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    [pathstr2, name2] = fileparts(pathstr);
    [pathstr3, name3] = fileparts(pathstr2);
    sFiles = {fullfile(name3, name2, [name, '.mat'])};
    idx = strfind(name,'ec'); idx1 = strfind(name,'_');
    sub_sel = name(idx+2:idx+5);
    run_sel = name(idx1(end-3)+1:idx1(end-2)-1);
    iSubject = find(contains(unq_bs_subj, sub_sel)==1);
    cd(pathstr2)
    if  ~isfolder(name(11:end))
        disp(sFiles)
%         pause,
        import_raw_to_db(sFiles{1});
    end
end

%% Est. head model
d1 = rdir(fullfile(BS_data_dir,'/*/ec*raw_low_clean/channel_vectorview306_acc1.mat'));
d2 = rdir(fullfile(BS_data_dir,'/*/ec*raw_tsss_low_clean/channel_vectorview306_acc1.mat'));
d3 = rdir(fullfile(BS_data_dir,'/*/EC*raw_low_clean/channel_vectorview306_acc1.mat'));
d4 = rdir(fullfile(BS_data_dir,'/*/ec*raw_clean_low/channel_vectorview306_acc1.mat'));
d5 = rdir(fullfile(BS_data_dir,'/*/ec*elecfix_raw_clean*/channel_vectorview306_acc1.mat'));

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

%% Clean-up (ICA + reject trials)
if analysis_flag.ica == 1
addpath(ft_path);
ft_defaults

d = rdir(fullfile(BS_data_dir,'/*/ec*raw_low_clean/channel_vectorview306_acc1.mat'));

for i=1:length(d)
    
    [pathstr, name] = fileparts(d(i).name);
    cd(pathstr)
    
    dd_results_3 = rdir(fullfile(pathstr,'results*.mat'));
    results_comments = [];
    for jj=1:length(dd_results_3)
        tmp = load(dd_results_3(jj).name);
        results_comments{jj} = tmp.Comment(1);
    end
    
    dd_32 = rdir(fullfile(pathstr,'data_3*_trl.mat'));
    if isempty(dd_32)
        dd_3 = rdir(fullfile(pathstr,'data_3*.mat'));
        sFiles = [];
        for j=1:length(dd_3)
            [pathstr, name] = fileparts(dd_3(j).name);
            [pathstr2, name2] = fileparts(pathstr);
            [pathstr3, name3] = fileparts(pathstr2);
            sFiles{j} = fullfile(name3, name2, [name, '.mat']);
        end
%         disp(sFiles(i))
        idx = strfind(name,'_');
        if ~isempty(results_comments)
            runidx = contains(results_comments, name(idx(1)+1));
        else
            runidx = [];
        end
        if isempty(runidx)
            aks_ica = input('run ICA (yes=1)?'); 
            if aks_ica ==1 
                bs_preprocess_ICA_reject_trials(sFiles),
            else
                disp('done')
            end
        end
    end
    
    dd_22 = rdir(fullfile(pathstr,'data_2*_trl.mat'));
    if isempty(dd_22)
        dd_2 = rdir(fullfile(pathstr,'data_2*.mat'));
        sFiles = [];
        for j=1:length(dd_2)
            [pathstr, name] = fileparts(dd_2(j).name);
            [pathstr2, name2] = fileparts(pathstr);
            [pathstr3, name3] = fileparts(pathstr2);
            sFiles{j} = fullfile(name3, name2, [name, '.mat']);
        end
%         disp(sFiles(i))
        idx = strfind(name,'_');
        if ~isempty(results_comments)
            runidx = contains(results_comments, name(idx(1)+1));
        else
            runidx = [];
        end
        if isempty(runidx)
            aks_ica = input('run ICA (yes=1)?');
            if aks_ica ==1 
                bs_preprocess_ICA_reject_trials(sFiles),
            else
                disp('done')
            end
        end
    end
end
end

%% Source modelling, DICS-BF
% rerun = {'EC1081','EC1082', 'EC1102', 'EC1121', 'EC1125', 'EC1154'};

d1 = rdir(fullfile(BS_data_dir,'/*/ec*raw_low_clean/channel_vectorview306_acc1.mat'));
d2 = rdir(fullfile(BS_data_dir,'/*/ec*raw_tsss_low_clean/channel_vectorview306_acc1.mat'));
d3 = rdir(fullfile(BS_data_dir,'/*/EC*raw_low_clean/channel_vectorview306_acc1.mat'));
d4 = rdir(fullfile(BS_data_dir,'/*/ec*raw_clean_low/channel_vectorview306_acc1.mat'));
d5 = rdir(fullfile(BS_data_dir,'/*/ec*elecfix_raw_clean*/channel_vectorview306_acc1.mat'));

d = [d1;d2;d3; d4; d5];

for i=1:length(d)
    
    [pathstr, name] = fileparts(d(i).name);
    
    dd_results = rdir(fullfile(pathstr,'results_subtraction*.mat'));
    results = [];
    for j=1:length(dd_results)
        cmt = load(dd_results(j).name);
        results{j} = cmt.Comment;
    end
    
    dd_32 = rdir(fullfile(pathstr,'data_3*_trl.mat'));
    if isempty(dd_32)
        dd_32 = rdir(fullfile(pathstr,'data_3_ic*.mat'));
        if isempty(dd_32)
            dd_32 = rdir(fullfile(pathstr,'data_3_trial*.mat'));
        end
    end
    
    dd_22 = rdir(fullfile(pathstr,'data_2*_trl.mat'));
    if isempty(dd_22)
        dd_22 = rdir(fullfile(pathstr,'data_2_ic*.mat'));
        if isempty(dd_22)
            dd_22 = rdir(fullfile(pathstr,'data_2_trial*.mat'));
        end
    end
    
    if isempty(results)
        flag = 1;
    elseif isempty(find(contains(results, 'wDICS')==1, 1))
        flag = 1;
    else
        flag = 0;
    end
    
%     if contains(pathstr, rerun)
%         flag = 1;
%     end
    
    if ~isempty(dd_32) && ~isempty(dd_32) && flag == 1 
        
        sFiles = [];
        for j=1:length(dd_32)
            [pathstr, name] = fileparts(dd_32(j).name);
            [pathstr2, name2] = fileparts(pathstr);
            [~, name3] = fileparts(pathstr2);
            sFiles{j} = fullfile(name3, name2, [name, '.mat']);
        end
        
        sFiles2 = [];
        for j=1:length(dd_22)
            [pathstr, name] = fileparts(dd_22(j).name);
            [pathstr2, name2] = fileparts(pathstr);
            [~, name3] = fileparts(pathstr2);
            sFiles2{j} = fullfile(name3, name2, [name, '.mat']);
        end
        
        disp(pathstr)
        
        % Process: FieldTrip: ft_sourceanalysis (wDICS) window contrast
        sFiles = bst_process('CallProcess', 'process_ft_sourceanalysis_dics_wcontrast', sFiles, sFiles2, ...
            'sensortype', 'MEG', ...  % MEG
            'poststim',   [7.21644966e-16, 2], ...
            'foi',        18, ...
            'tpr',        4, ...
            'tlength',    0.3, ...
            'ovp',        0.5, ...
            'simaps',     1, ...
            'avmaps',     1, ...
            'method',     'subtraction', ...  % Subtraction (post-pre)
            'erds',       'erd', ...  % ERD
            'effect',     'abs', ...  % abs
            'maxfreq',    40, ...
            'showtfr',    1);
    end
end

%% Intra-subject averaging
load(protocol);
Subj = ProtocolSubjects.Subject;
ProtocolSubjects = []; k=1;
for i=1:length(Subj)
    if contains(Subj(i).Name, 'EC')
        ProtocolSubjects{k} =  Subj(i).Name; k=1+k;
    end
end

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
    'EC1111'};

incomplete_data = {'EC1127'}; % only one run
subj = unq_bs_subj;

clc
close all,
subj_del = [];
for ii = 1:length(subj)
    
    if ~(contains(subj{ii},no_anat)) || contains(subj{ii},incomplete_data)
        cd(fullfile(BS_data_dir,subj{ii}))
        dd = rdir(['./*',tag,'*/results*.mat']);
        for jj=1:length(dd), disp([num2str(jj),':',dd(jj).name]); end
        sel = [];
        for jj=1:length(dd)
            tmp = load(dd(jj).name);
            disp(tmp.Comment);
            if contains(tmp.Comment,'wDICS: subtraction')
                sel = [sel,jj];
            end
        end
        
        dd_Sel = [];
        if length(sel)>2
            clc
            subj_del = [subj_del;subj{ii}];
            dd_Sel = dd(sel);
            for jj=1:length(dd_Sel)
                tmp = load(dd_Sel(jj).name);
                disp([num2str(jj),':',dd_Sel(jj).name])
                disp(tmp.Comment);
            end
            delin = input('del extra files:');
            delete(dd_Sel(delin).name)
            %
            dd = rdir(['./*',tag,'*/results*.mat']);
            for jj=1:length(dd), disp([num2str(jj),':',dd(jj).name]); end
            sel = [];
            for jj=1:length(dd)
                tmp = load(dd(jj).name);
                disp(tmp.Comment);
                if contains(tmp.Comment,'s_2sided')
                    sel = [sel,jj];
                end
            end
        end
        
        dd_Sel = [];
        if length(sel)==2
            dd_Sel = dd(sel);
            
            sFiles1 = []; sFiles_name = [];
            for jj=1:length(dd_Sel)
                sFiles1{jj} = fullfile(subj{ii},dd_Sel(jj).name);
                [aa,bb] = fileparts(dd(jj).name);
                sFiles_name{jj} = aa;
            end
            disp(sFiles1')
            disp('------------');
            
            dd1 = rdir('./@intra/results*.mat');
            if ~isempty(dd1)
                for kk=1:length(dd1)
                    tmp = load(dd1(kk).name);
                    disp(tmp.Comment);
                    if contains(tmp.Comment,[subj{ii}, '_wDICS_'])
                        runok = 0; break,
                    else
                        runok = 1;
                    end
                end
            else
                runok = 1;
            end
            if (length(sFiles1)==2
                % Process: Average: Everything
                bst_process('CallProcess', 'process_average', sFiles1, [], ...
                    'avgtype',         1, ...  % Everything
                    'avg_func',        1, ...  % Arithmetic average:  mean(x)
                    'weighted',        0, ...
                    'Comment', [subj{ii}, '_wDICS_contrast'], ...
                    'scalenormalized', 0);
            else
                warning(['check data:', subj{ii}])
            end
        end
    end
end
disp('intra-subject source averaging was completed!');

%%
cd(BS_data_dir);
dd = rdir(fullfile (BS_data_dir,'/*/@intra/results*.mat'));
dd1 = rdir(fullfile (BS_data_dir,'/Group*/@intra/results*.mat'));

d = []; d2 = [];
for ii=1:length(dd), d{ii}=dd(ii).name; disp(dd(ii).name); end
for ii=1:length(dd1), d2{ii}=dd1(ii).name; disp(dd1(ii).name); end
dd3 = setdiff(d,d2);

% atag = {'3','2'};

subj = []; comm_data = []; sel = []; dconn = []; kk = 1;
sFiles_name = [];
for ii=1:length(dd3)
    [pp, ~] = fileparts(dd3{ii});
    disp([num2str(ii),'/',num2str(length(dd3))])
    cd(pp)
    tmp = load(dd3{ii});
    comm_data{ii} = tmp.Comment;
    if contains(comm_data{ii}, ['_wDICS_contrast'])
        %             pause,
        sel = [sel,ii];
        tkz = tokenize(comm_data{ii},'_');
        subj{kk}= tkz{1}(6:end);
        dconn{kk} = tkz{2};
        %             pause
        sFiles_name{kk} = [subj{kk},'_',dconn{kk}];
        kk=1+kk;
    end
end
d_sel = dd3(sel);

% looking for already caculated maps ..
dd = rdir(fullfile (BS_data_dir,'/Group*/wDICS_contrast_18_4/results*.mat'));
subj_comp = []; comm_data = []; sel = []; dconn = []; kk = 1;
sFiles_name_completed = [];
for ii=1:length(dd)
    [pp, ~] = fileparts(dd(ii).name);
    disp([num2str(ii),'/',num2str(length(dd))])
    cd(pp)
    tmp = load(dd(ii).name);
    comm_data{ii} = tmp.Comment;
    disp(comm_data{ii})
    if contains(comm_data{ii}, ['_wDICS_contrast'])
%         pause,
        sel = [sel,ii];
        tkz = tokenize(comm_data{ii},'_');
        subj_comp{kk}= tkz{1}(13:end);
        dconn{kk} = tkz{2};
        sFiles_name_completed{kk} = [subj_comp{kk}];
        kk=1+kk;
    end
end
if isempty(sFiles_name_completed)
    sFiles_name_completed = 'none';
end

% Project on default anatomy (for group mapping)
idx = find(contains(sFiles_name,'wDICS')==1);
destSurfFile = '@default_subject/tess_cortex_pial_low.mat';

load(protocol);
Subj = ProtocolSubjects.Subject;
ProtocolSubjects = []; k=1;
for i=1:length(Subj)
    if contains(Subj(i).Name, 'EC')
        ProtocolSubjects{k} =  Subj(i).Name; k=1+k;
    end
end

cd(BS_data_dir);
sFiles = [];
if length(idx)>1
    for i=1:length(idx)
        %             pause,
        if ~contains(subj{i},no_anat) && ~contains(sFiles_name{i},sFiles_name_completed)
            disp(i);
            sFiles = d_sel{i};
            tkz = tokenize(sFiles,'/');
            idx1 = find(strcmp(tkz{end-2}, ProtocolSubjects)==1);
            db_reload_conditions(idx1);
            db_reload_subjects(idx1);
            disp(tkz{end-2});
            sFiles1 = fullfile(tkz{end-2}, tkz{end-1}, tkz{end});
            pause,
            bst_project_sources ({sFiles1}, destSurfFile);
        end
    end
end



