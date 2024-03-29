%% Logopenicppa dataset, Medical College of Wisconsin

% Script: BS Process (preprocessing, source analysis)
% Project: Logopenicppa_CRM
% Writtern by: Vahab YoussofZadeh
% Update: 01/04/2023

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
cd('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Logopenicppa')

indir = '/data/MEG/Research/logopenicppa/raw';
outdir = '/data/MEG/Research/logopenicppa/BS_process';

ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
addpath(ft_path); ft_defaults

%%
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/brewermap')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/Colormaps-from-MatPlotLib2.0')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/Miscellaneous')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Logopenicppa/functions')

%%
% adding BS path
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3'; %BS-2021
addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

BS_dir = '/data/MEG/Research/logopenicppa/BS_process/Logopenicppa_CRM_new';
BS_data_dir = fullfile(BS_dir,'data');
protocol = fullfile(BS_dir, 'data/protocol.mat');

%% Loading up raw data
tag = 'CRM'; %- Semantic decision

cfg = [];
cfg.indir = indir;
cfg.tag = tag;
datafile = do_datalookup(cfg);

%%
% Input files
sFiles = {...
    'link|32037_001_1/Run04_CRM_run2_raw_clean_low/results_PNAI_MEG_GRAD_MEG_MAG_KERNEL_220309_1220.mat|32037_001_1/Run04_CRM_run2_raw_clean_low/data_1_average_220309_1220.mat', ...
    'link|32037_001_1/Run03_CRM_run1_raw_clean_low/results_PNAI_MEG_GRAD_MEG_MAG_KERNEL_220309_1431.mat|32037_001_1/Run03_CRM_run1_raw_clean_low/data_1_average_220309_1431.mat'};

% Start a new report
bst_report('Start', sFiles);

% Process: Average: Everything
sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'scalenormalized', 0);

%% intra sub avg.
destSurfFile = '@default_subject/tess_cortex_pial_low.mat';

load(protocol);
Subj = ProtocolSubjects.Subject;
ProtocolSubjects = []; k=1;
for i=1:length(Subj)
    if contains(Subj(i).Name, '32037')
        ProtocolSubjects{k} =  Subj(i).Name; k=1+k;
    end
end

subj = ProtocolSubjects;
tag = 'CRM*raw_clean_low'; %- Semantic decision

datatag = {'1'};
for ii = 1:length(subj)
    cd(fullfile(BS_data_dir,subj{ii}))
    
    dd1 = rdir('./@intra/results_average_*.mat');
    Comments = [];
    for jj=1:length(dd1)
        tmp = load(dd1(jj).name);
        Comments{jj} = tmp.Comment;
    end
    
    if isempty(find(contains(Comments, subj{ii})==1, 1))        
        dd = rdir(['./*',tag,'*/results_PNAI*.mat']);
        for jj=1:length(dd), disp([num2str(jj),':',dd(jj).name]); end
        sFiles1 = []; sFiles_name = [];
        for jj=1:length(dd)
            sFiles1{jj} = fullfile(subj{ii},dd(jj).name);
            [aa,bb] = fileparts(dd(jj).name);
            sFiles_name{jj} = aa;
        end
        
        sFiles2 = [];
        k=1;
        for j=1:length(sFiles1)
            cd(BS_data_dir)
            tmp  = load(sFiles1{j});
            disp(tmp.Comment)
            d = rdir(['./', fileparts(sFiles1{j}), '/data_', datatag{1}, '_average*.mat']);
            disp(sFiles1{j})
            if length(d) > 1
                for jj=1:length(d)
                    tmp2 = load(d(jj).name);
                    disp([num2str(jj), ': ', tmp2.Comment])
                end
                del_sel = input('sel_del:');
                delete(d(del_sel).name)
                d = rdir(['./',fileparts(sFiles1{j}), '/data_', datatag{dtagsel}, '_average*.mat']);
            end
            if ~contains(tmp.Comment, 'ms')
                sFiles2{k} = ['link|',sFiles1{j}, '|',d.name(3:end)]; k=k+1;
            end
        end
        
        disp(sFiles2')
        disp(subj{ii})
        % Process: Average: Everything
        bst_process('CallProcess', 'process_average', sFiles2, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...
            'Comment',        ['Avg:', subj{ii}], ...
            'scalenormalized', 0);
    end
end


%%
cd(BS_data_dir)
dd = rdir(fullfile('./Group_analysis/*/results_*.mat'));

Comments_all = [];
for jj=1:length(dd)
    tmp = load(dd(jj).name);
    Comments_all{jj} = tmp.Comment;
end


for ii = 1:length(subj)
    disp(subj{ii})
    cd(fullfile(BS_data_dir,subj{ii}))    
    dd1 = rdir('./@intra/results_average_*.mat');
    Comments = [];
    for jj=1:length(dd1)
        tmp = load(dd1(jj).name);
        Comments{jj} = tmp.Comment;
    end
    
    if ~contains(Comments_all, subj{ii})
        idx = contains(Comments, subj{ii})==1;
%         pause,
        sFiles = bst_project_sources({fullfile(subj{ii}, dd1(idx).name)}, destSurfFile, 0, 1);
    end
end

%% Loading up raw data
% cd(indir)
% 
% clc
% if exist(['datalog_',tag,'.mat'], 'file') == 2
%     load(['datalog_',tag,'.mat'])
% else
%     clear datafolder datafile subj_all sub
%     datafile_fif = [];
%     d1 = rdir([indir,['/**/tsss/','*CRM*_raw.fif']]);
% %     d2 = rdir([indir,['/**/tsss/*',tag,'*_raw_tsss.fif']]);
%     d = d1;
%     for i=1:length(d)
%         [pathstr, name] = fileparts(d(i).name);
%         datafolder{i} = pathstr;
%         datafile{i} = d(i).name;
%         Index = strfind(datafile{i}, '/');
%         sub{i} = datafile{i}(Index(6)+1:Index(7)-1);
%         subj_all{i} = [num2str(i), ': ', datafile{i}(Index(6)+1:Index(7)-1)];
%     end
%     datafile_fif = vertcat(datafile_fif,datafile);
%     datafile_fif = datafile_fif';
%     save(['datalog_',tag,'.mat'],'datafile_fif','subj_all', 'sub')
% end
% disp(datafile_fif)

%%
% d = rdir([BS_data_dir,'/*/Run*_low/results*KERNEL*.mat']);
% 
% clear subj subj_data unq_bs_subj
% for i=1:length(d)
%     %     tmp = load(d(i).name);
%     tkz  = tokenize(d(i).name,'/');
%     subj{i} = tkz{end-2};
%     subj_data{i} = tkz{end-1};
% end
% unq_bs_subj = unique(subj);

%%
% destSurfFile = '@default_subject/tess_cortex_pial_low.mat';
% 
% clc
% clear nsFiles
% for  i=1:length(unq_bs_subj)
%     idx = find(contains(subj, unq_bs_subj{i})==1);
%     if length(idx)==2
%         for j=1:length(idx)
%             D_name = d(idx(j)).name;
%             tmp = load(D_name);
%             disp(D_name);
%             idx2 = strfind(D_name,'/data/');
%             [pathstr1, name1] = fileparts(D_name);
%             d2 = rdir([pathstr1,'/data_1_average*.mat']); 
%             [pathstr2, name2] = fileparts(d2.name);
%             [pathstr3, name3] = fileparts(D_name(idx2(2)+6:end));
%             sFiles{j} = ['link|', D_name(idx2(2)+6:end),'|', pathstr3, '/', name2,'.mat'];
%             nsFiles(j) = bst_project_sources({sFiles{j}}, destSurfFile, 0, 1);
%         end
%         % Process: Average: Everything
%         bst_process('CallProcess', 'process_average', nsFiles, [], ...
%             'avgtype',         1, ...  % Everything
%             'avg_func',        1, ...  % Arithmetic average:  mean(x)
%             'weighted',        0, ...
%             'Comment', unq_bs_subj{i}, ...
%             'scalenormalized', 0);
% %         pause,
%     end
% end

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
no_anat = {''};
sub_all1 = setdiff(unq_bs_subj,no_anat);

%% Import raw data
% L_data =length(datafile_fif);
% not_imported = [];
% for ii=1:L_data
%     [~, name] = fileparts(datafile_fif{ii});
%     idx = strfind(name,'_'); 
%     sub_sel = name(3:idx(1)-1); 
%     run_sel = name(idx(2)+1:idx(3)-1);
%     
%     datadir_sub = fullfile(BS_dir,'data/',['32037_',sub_sel]);
%     cd(datadir_sub)
%     idx = strfind(datadir_sub,'/');
%     okrun = find(contains(no_anat,datadir_sub(idx(end)+1:end))==1);
%     if  ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw']) ...
%             && isempty(okrun) && ...
%             ~isfolder(['@rawEC',sub_sel, '_SD_', run_sel, '_raw'])...
%             && ~isfolder(['@rawEC',sub_sel, '_SD_', run_sel, '_elecfix_raw']) ...
%             && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_elecfix_raw'])...
%             && ~isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw_tsss'])
%         iSubject = find(contains(unq_bs_subj, sub_sel)==1);
%         RawFiles = datafile_fif{ii};
%         pause(3),
%         OutputFiles = import_raw(RawFiles, 'FIF', iSubject, [], []);
%     end
% end

%% Preprocess raw data
% d1 = rdir(fullfile(BS_data_dir,'/EC*/@raw*/*_raw.mat'));
% d2 = rdir(fullfile(BS_data_dir,'/EC*/@raw*/*_raw_tsss.mat'));
% d = [d1;d2];
% for i=1:length(d)
%     [pathstr, name] = fileparts(d(i).name);
%     [pathstr2, name2] = fileparts(pathstr);
%     [pathstr3, name3] = fileparts(pathstr2);
%     sFiles = {fullfile(name3, name2, [name, '.mat'])};
%     idx = strfind(name,'ec');
%     if isempty(idx)
%         idx = strfind(name,'EC');
%     end
% %     idx1 = strfind(name,'_');
%     idx2 = strfind(name,'run');
%     if isempty(idx2)
%        idx2 = strfind(name,'Run'); 
%     end
%     sub_sel = name(idx+2:idx+5);
%     run_sel = name(idx2+3);
%     cd(pathstr2)
%     if  ~isfolder(['@rawec',sub_sel, '_SD_run', run_sel, '_raw_low']) ...
%             && ~isfolder(['@rawEC',sub_sel, '_SD_run', run_sel, '_raw_low'])...
%             && ~isfolder(['@rawEC',sub_sel, '_SD_run', run_sel, '_elecfix_raw_low'])...
%             && ~isfolder(['@rawec',sub_sel, '_SD_run', run_sel, '_elecfix_raw_low'])...
%             && ~isfolder(['@rawec',sub_sel, '_SD_Run', run_sel, '_raw_low'])...
%             && ~isfolder(['@rawec',sub_sel, '_SD_Run', run_sel, '_raw_tsss_low'])...
%             && ~isfolder(['@rawec',sub_sel, '_SD_run', run_sel, '_raw_tsss_low'])
%         disp(sFiles)
% %         pause,
% %         cd(BS_data_dir)
%         bs_preprocess(sFiles)
%     end
% end

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

 %%
% d1 = rdir(fullfile(BS_data_dir,'/**/@raw*/*raw_low_clean.mat'));
% d2 = rdir(fullfile(BS_data_dir,'/**/@raw*/*raw_tsss_low_clean.mat'));
% d = [d1;d2];
% for i=1:length(d)
%     [pathstr, name] = fileparts(d(i).name);
%     [pathstr2, name2] = fileparts(pathstr);
%     [pathstr3, name3] = fileparts(pathstr2);
%     sFiles = {fullfile(name3, name2, [name, '.mat'])};
%     idx = strfind(name,'ec'); idx1 = strfind(name,'_');
%     sub_sel = name(idx+2:idx+5);
%     run_sel = name(idx1(end-3)+1:idx1(end-2)-1);
%     iSubject = find(contains(unq_bs_subj, sub_sel)==1);
%     cd(pathstr2)
%     if  ~isfolder(name(11:end))
%         disp(sFiles)
% %         pause,
%         import_raw_to_db(sFiles{1});
%     else
%         dd = rdir(fullfile(name3, name2, name(11:end), 'data_*_average*.mat'));
%         cd(fullfile(name3, name2, name(11:end)))
%     end
%     
% end

%% Est. head model
d1 = rdir(fullfile(BS_data_dir,'/*/ec*raw_low_clean/channel_vectorview306_acc1.mat'));
d2 = rdir(fullfile(BS_data_dir,'/*/ec*raw_tsss_low_clean/channel_vectorview306_acc1.mat'));
d3 = rdir(fullfile(BS_data_dir,'/*/EC*raw_low_clean/channel_vectorview306_acc1.mat'));
d = [d1;d2;d3];

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

%% Source, lcmv
d1 = rdir(fullfile(BS_data_dir,'/*/ec*raw_low_clean/channel_vectorview306_acc1.mat'));
d2 = rdir(fullfile(BS_data_dir,'/*/ec*raw_tsss_low_clean/channel_vectorview306_acc1.mat'));
d3 = rdir(fullfile(BS_data_dir,'/*/EC*raw_low_clean/channel_vectorview306_acc1.mat'));
d = [d1;d2;d3];

for i=1:length(d)
    
    [pathstr, name] = fileparts(d(i).name);
    dd_results = rdir(fullfile(pathstr,'results*.mat'));
    results = [];
    for j=1:length(dd_results)
        cmt = load(dd_results(j).name);
        results{j} = cmt.Comment(1);
    end
    
    dd_32 = rdir(fullfile(pathstr,'data_3*_trl.mat'));
    if isempty(dd_32)
        dd_32 = rdir(fullfile(pathstr,'data_3_ic*.mat'));
        if isempty(dd_32)
            dd_32 = rdir(fullfile(pathstr,'data_3_trial*.mat'));
        end
    end
    
    if isempty(results)
        flag = 1;
    elseif isempty(find(contains(results, '3')==1, 1))
        flag = 1;
    else
        flag = 0;
    end
    
    if ~isempty(dd_32) && flag == 1
        sFiles = [];
        for j=1:length(dd_32)
            [pathstr, name] = fileparts(dd_32(j).name);
            [pathstr2, name2] = fileparts(pathstr);
            [~, name3] = fileparts(pathstr2);
            sFiles{j} = fullfile(name3, name2, [name, '.mat']);
        end
        cfg = []; cfg.comment = '3'; bs_lcmv(cfg,sFiles)
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
    elseif isempty(find(contains(results, '2')==1, 1))
        flag = 1;
    else
        flag = 0;
    end
    
    if ~isempty(dd_22) && flag == 1
        sFiles = [];
        for j=1:length(dd_22)
            [pathstr, name] = fileparts(dd_22(j).name);
            [pathstr2, name2] = fileparts(pathstr);
            [~, name3] = fileparts(pathstr2);
            sFiles{j} = fullfile(name3, name2, [name, '.mat']);
        end
        cfg = []; cfg.comment = '2'; bs_lcmv(cfg, sFiles)
    end
end

%%
% d = rdir(fullfile(BS_data_dir,'/*/ec*raw_low_clean/channel_vectorview306_acc1.mat'));
% 
% clear datafolder datafile subj_all run_all
% for i=1:length(d)
%     [pathstr, name] = fileparts(d(i).name);
%     datafolder{i} = pathstr;
%     datafile{i} = d(i).name;
%     Index = strfind(datafile{i}, 'EC'); Index1 = strfind(datafile{i}, '/');
%     Index2 = strfind(datafile{i}, '_');
%     subj_all{i} = datafile{i}(Index(2):Index1(end-1)-1);
%     run_all{i} = datafile{i}(Index2(end-5)+1:Index2(end-4)-1);
% end
% BS_chan = [];
% BS_chan.run_all = run_all;
% BS_chan.subj_all = subj_all;
% BS_chan.datafile = datafile;
% BS_chan.datafolder = datafolder;
% BS_datafile = datafile;

%% Project to default template
destSurfFile = '@default_subject/tess_cortex_pial_low.mat';

load(protocol);
Subj = ProtocolSubjects.Subject;
ProtocolSubjects = []; k=1;
for i=1:length(Subj)
    if contains(Subj(i).Name, 'EC')
        ProtocolSubjects{k} =  Subj(i).Name; k=1+k;
    end
end

subj = ProtocolSubjects;

datatag = {'3','2'};
for ii = 1:length(subj)
    cd(fullfile(BS_data_dir,subj{ii}))
    dd = rdir(['./*',tag,'*/results_PNAI*.mat']);
    for jj=1:length(dd), disp([num2str(jj),':',dd(jj).name]); end
    sFiles1 = []; sFiles_name = [];
    for jj=1:length(dd)
        sFiles1{jj} = fullfile(subj{ii},dd(jj).name);
        [aa,bb] = fileparts(dd(jj).name);
        sFiles_name{jj} = aa;
    end
    for j=1:length(sFiles1)
        cd(BS_data_dir)
        tmp  = load(sFiles1{j}); disp(tmp.Comment(1))
        dtagsel = find(contains(datatag,tmp.Comment(1))==1);
        d = rdir(['./',fileparts(sFiles1{j}), '/data_', datatag{dtagsel}, '_average*.mat']);
        disp(sFiles1{j})
        if length(d) > 1
            for jj=1:length(d)
                tmp2 = load(d(jj).name);
                disp([num2str(jj), ': ', tmp2.Comment])
            end
            del_sel = input('sel_del:');
            delete(d(del_sel).name)
            d = rdir(['./',fileparts(sFiles1{j}), '/data_', datatag{dtagsel}, '_average*.mat']);
        end
        sFiles2 = ['link|',sFiles1{j}, '|',d.name(3:end)];
        
        tkz = tokenize(sFiles1{j},'/');
        d_name = [];
        dd = rdir(fullfile(BS_data_dir, '/Group_analysis', tkz{2}, 'results_*.mat'));
        
        for jj = 1:length(dd)
            d_name{jj} = dd(jj).name;
        end
        tkz2 = tokenize(sFiles1{j},'_');
        
        % --- Condition OK
        nn = d.name(3:end);
        idd = strfind(nn,'data_');
        
        cond_ok = [];
        if ~isempty(d_name)
            tmp_comment = [];
            for k=1:length(d_name)
                tmp  = load(d_name{k});
                tmp_comment{k} = tmp.Comment;
            end
            cond_ok = double(isempty(find(contains(tmp_comment,[' ',nn(idd+5),' ']), 1)));
        else
            cond_ok = 0;
        end
        %---
        
        if isempty(d_name)
            run_ok =1;
        elseif ~contains(d_name, tkz2{end}(1:end-4))
            run_ok =1;
        else
            run_ok = 0;
        end
        
        if (run_ok ==1) || (cond_ok == 1)
            pause
%             sFiles = bst_process('CallProcess', 'process_project_sources', sFiles2, [], ...
%                 'headmodeltype', 'surface');  % Cortex surface
            sFiles = bst_project_sources({sFiles2}, destSurfFile, 0, 1);
        end
    end
end

%% delete extra smoothed files
cd(BS_data_dir)
clc
dd = rdir(fullfile('./Group_analysis/*/results_*.mat'));

for ii=1:length(dd)
    disp([num2str(ii), '/', num2str(length(dd))])
    cd(BS_data_dir)
    [path, name] = fileparts(dd(ii).name);
    cd(path)
    d = rdir('./*_ssmooth*.mat');
    
    comment = [];
    for j=1:length(d)
        tmp = load(d(j).name);
%         disp([num2str(j), ': ', tmp.Comment])
        comment{j} = tmp.Comment(21);
    end
    
    if length(d) > 2 || (length(d) == 2 && comment{1} == comment{2})
        for j=1:length(d)
            tmp = load(d(j).name);
            disp([num2str(j), ': ', tmp.Comment])
        end
        del_sel = input('sel to delete:');
        delete(d(del_sel).name)
    end
end

%% delete extra avg file
cd(BS_data_dir)
clc
dd = rdir(fullfile('./Group_analysis/*/results_*.mat'));

for ii=1:length(dd)
    disp([num2str(ii), '/', num2str(length(dd))])
    cd(BS_data_dir)
    [path, name] = fileparts(dd(ii).name);
    cd(path)
    
    if ~contains(path,'@intra')
%         pause,
        d = rdir('./*results*.mat');
        
        comment = []; k=1;
        for j=1:length(d)
         tmp = load(d(j).name);
            if ~contains(tmp.Comment,'ssmooth') && ~contains(tmp.Comment,'Avg: ')
                comment{k} = tmp.Comment(13); k=1+k;
            end
        end
%         if tmp.Time(end) <= 0.30
%             
%         end
        
        if length(comment) > 2 
            for j=1:length(d)
                tmp = load(d(j).name);
                disp([num2str(j), ': ', tmp.Comment])
            end
            del_sel = input('sel to delete:');
            delete(d(del_sel).name)
        end
    end
end

%% delete short epochs
cd(BS_data_dir)
clc
dd = rdir(fullfile('./EC*/*/results_*.mat'));

for ii=1:length(dd)
    
    cd(BS_data_dir)
    [path, name] = fileparts(dd(ii).name);
    cd(path)
    
    if ~contains(path,'@intra')
        d = rdir('./*_average*.mat');
        comment = []; tt=[];
        for j=1:length(d)
            tmp = load(d(j).name);
            tt(j) = tmp.Time(end);
            if round(tt(j)) ~= 2
                rmv = 1; break,
            else
                rmv = 0;
            end
        end
        
        if rmv ==1
            disp(path)
%             pause
%             rmdir(fullfile(BS_data_dir, path),'s')
        end
    end
end

%% Apply smoothing, for group source analysis
cd(BS_data_dir)
clc
dd = rdir(fullfile('./Group_analysis/*clean/results_*.mat'));
ft_progress('init', 'text',     'please wait ...');
comment = [];
for j=1:length(dd)
    ft_progress(i/length(dd), 'Processing event %d from %d', j, length(dd));
    tmp = load(dd(j).name); comment{j} = tmp.Comment;
    pause(0.1);
end
ft_progress('close');

%%
% ft_progress('init', 'text',     'please wait ...');
comm_data = [];
for ii=1:length(dd) 
    cd(BS_data_dir)
    disp([num2str(ii), '/', num2str(length(dd))])
%     ft_progress(i/length(dd), 'Processing event %d from %d', j, length(dd));
    [a, ~] = fileparts(dd(ii).name);
    tmp = load(dd(ii).name);
    comm_data_sel = tmp.Comment;
    if contains(comm_data_sel,datatag) && ~contains(comm_data_sel,'ssmooth') ...
            && ~contains(comm_data_sel,'matlab') && ...
            ~contains(comm_data_sel,' - ') && ...
            ~contains(comm_data_sel,'mean_')
        disp(ii)
        disp(comm_data_sel);
        tkz = tokenize(dd(ii).name,'/');
        sFiles = dd(ii).name(3:end);
        
        [path, name] = fileparts(dd(ii).name);
        ddd = rdir([path, '/results_PNAI*.mat']);
        comment1 = [];
        for j=1:length(ddd)
            tmp = load(ddd(j).name);
            comment1{j} = tmp.Comment;
        end
        if ~isempty(comment1)
            if isempty(find(contains(comment1, ['ssmooth', '_', comm_data_sel])==1, 1))
                %             disp(comm_data_sel)
                %             [path, name] = fileparts(dd(ii).name);
                %             cd(path)
                %             d = rdir('./*_ssmooth*.mat');
                %             label = [];
                %             if length(d) > 1
                %                 for j=1:length(d)
                %                     tmp = load(d(j).name);
                %                     disp([num2str(j), ': ', tmp.Comment])
                %                     label{j} = tmp.Comment(21);
                %                 end
                %                 if label{1} == label{2}
                %                     del_sel = input('sel to delete:');
                %                     delete(d(del_sel).name)
                %                 end
                %             end
                %             pause,
                bst_process('CallProcess', 'process_ssmooth_surfstat', sFiles, [], ...
                    'fwhm',       3, ...
                    'overwrite',  0, ...
                    'Comment', ['ssmooth', '_', comm_data_sel], ...
                    'source_abs', 1);
            end
        end
    end
end
ft_progress('close');


%% Intra-subject averaging
% db_reload_database('current',1)
% load(protocol);
% Subj = ProtocolSubjects.Subject;
% ProtocolSubjects = []; k=1;
% for i=1:length(Subj)
%     if contains(Subj(i).Name, 'EC')
%         ProtocolSubjects{k} =  Subj(i).Name; k=1+k;
%     end
% end

clc
close all,
subj_del = [];

cd(BS_data_dir)
dd = rdir(fullfile('./Group_analysis/*/results_*.mat'));
for jj=1:length(dd), disp([num2str(jj),':',dd(jj).name]); end

sFiles_name = [];
for jj=1:length(dd)
    sFiles_name{jj} = fullfile(dd(jj).name(3:end));
end

Comment = [];
for jj=1:length(sFiles_name)
    cd(BS_data_dir)
    tmp  = load(sFiles_name{jj});
    Comment{jj} = tmp.Comment;
end

idx_31 = find(contains(Comment, '3_')==1); idx_32 = find(contains(Comment, '3 (')==1);
idx_21 = find(contains(Comment, '2_')==1); idx_22 = find(contains(Comment, '2 (')==1);

idx3 = [idx_31, idx_32];
idx2 = [idx_21, idx_22];

idx_ssmoth = find(contains(Comment, 'ssmooth')==1);
idx_3_smooth = intersect(idx3,idx_ssmoth);
idx_2_smooth = intersect(idx2,idx_ssmoth);

% Process: Average: Everything
bst_process('CallProcess', 'process_average', sFiles_name(idx_3_smooth), [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment', 'mean_3', ...
    'scalenormalized', 0);

% Process: Average: Everything
bst_process('CallProcess', 'process_average', sFiles_name(idx_2_smooth), [], ...
    'avgtype',         1, ...  % Everything
    'avg_func',        1, ...  % Arithmetic average:  mean(x)
    'weighted',        0, ...
    'Comment', 'mean_2', ...
    'scalenormalized', 0);

%% Subject avg
sFiles_3 = sFiles_name(idx_3_smooth);
sFiles_2 = sFiles_name(idx_2_smooth);

L = length(sFiles_3); k = 1;
clear subjs_3
for i=1:length(sFiles_3)
%     disp([num2str(i), '/' , num2str(length(sFiles_3))])
    tmp = load(sFiles_3{i});
    subjs_3{k} = (tmp.Comment(9:14));
    disp([subjs_3{k}, ': ', num2str(i), '/' , num2str(length(sFiles_3))])
    k=1+k;
end
unq_bs_subj_3 = unique(subjs_3);

L = length(sFiles_2); k = 1;
clear subjs_2
for i=1:length(sFiles_2)
%     disp([num2str(i), '/' , num2str(length(sFiles_3))])
    tmp = load(sFiles_2{i});
    subjs_2{k} = (tmp.Comment(9:14));
    disp([subjs_2{k}, ': ', num2str(i), '/' , num2str(length(sFiles_2))])
    k=1+k;
end
unq_bs_subj_2 = unique(subjs_2);

%%
unq_bs_subj_3_sel = {'EC1092',  'EC1112'};
unq_bs_subj_2_sel = {'EC1092',  'EC1112'};

s_in = sFiles_name(idx_3_smooth);
for j=1:length(unq_bs_subj_3_sel)
    disp([num2str(j), '/', num2str(length(unq_bs_subj_3_sel))])
    idx = contains(subjs_3, unq_bs_subj_3_sel{j})==1;
    bst_process('CallProcess', 'process_average', s_in(idx), [], ...
        'avgtype',         1, ...  % Everything
        'avg_func',        1, ...  % Arithmetic average:  mean(x)
        'weighted',        0, ...
        'Comment', ['Anim_', unq_bs_subj_3_sel{j}], ...
        'scalenormalized', 0);
end

s_in = sFiles_name(idx_2_smooth);
for j=1:length(unq_bs_subj_2_sel)
    disp([num2str(j), '/', num2str(length(unq_bs_subj_2_sel))])
    idx = contains(subjs_2, unq_bs_subj_2_sel{j})==1;
    bst_process('CallProcess', 'process_average', s_in(idx), [], ...
        'avgtype',         1, ...  % Everything
        'avg_func',        1, ...  % Arithmetic average:  mean(x)
        'weighted',        0, ...
        'Comment', ['Symbol_', unq_bs_subj_2_sel{j}], ...
        'scalenormalized', 0);
end

%%

s_in = sFiles_name(idx_3_smooth);
for j=1:length(unq_bs_subj_3)
    disp([num2str(j), '/', num2str(length(unq_bs_subj_3))])
    idx = contains(subjs_3, unq_bs_subj_3{j})==1;
    bst_process('CallProcess', 'process_average', s_in(idx), [], ...
        'avgtype',         1, ...  % Everything
        'avg_func',        1, ...  % Arithmetic average:  mean(x)
        'weighted',        0, ...
        'Comment', ['Anim_', unq_bs_subj_3{j}], ...
        'scalenormalized', 0);
end

s_in = sFiles_name(idx_2_smooth);
for j=1:length(unq_bs_subj_2)
    disp([num2str(j), '/', num2str(length(unq_bs_subj_2))])
    idx = contains(subjs_2, unq_bs_subj_2{j})==1;
    bst_process('CallProcess', 'process_average', s_in(idx), [], ...
        'avgtype',         1, ...  % Everything
        'avg_func',        1, ...  % Arithmetic average:  mean(x)
        'weighted',        0, ...
        'Comment', ['Symbol_', unq_bs_subj_2{j}], ...
        'scalenormalized', 0);
end

%%
sub_3 = [];
d_in = sFiles_name(idx_3_smooth);
for i=1:length(d_in)
   sub_3{i} = d_in{i}(16:21); 
end
usub_3 = unique(sub_3)';


sub_2 = [];
d_in = sFiles_name(idx_2_smooth);
for i=1:length(d_in)
   sub_2{i} = d_in{i}(16:21); 
end
usub_2 = unique(sub_2)';

A = usub_2;
B = usub_3;
setdiff(lower(A), lower(B))

for i=1:length(sub_3)
    if length(find(contains(sub_3, sub_3(i))))~=2
        disp(sub_3(i))
    end
end

for i=1:length(sub_2)
    if length(find(contains(sub_2, sub_2(i))))~=2
        disp(sub_2(i))
    end
end

%%
sFiles = bst_process('CallProcess', 'process_test_parametric1',  sFiles_name(idx_3_smooth), [], ...
    'timewindow',    [-0.3, 2], ...
    'scoutsel',      {}, ...
    'scoutfunc',     1, ...  % Mean
    'isnorm',        0, ...
    'avgtime',       0, ...
    'Comment',       'stat_3', ...
    'test_type',     'ttest_onesample', ...  % One-sample Student's t-test    X~N(m,s)t = mean(X) ./ std(X) .* sqrt(n)      df=n-1
    'tail',          'one+');  % One-tailed (+)

sFiles = bst_process('CallProcess', 'process_test_parametric1',  sFiles_name(idx_2_smooth), [], ...
    'timewindow',    [-0.3, 2], ...
    'scoutsel',      {}, ...
    'scoutfunc',     1, ...  % Mean
    'isnorm',        0, ...
    'avgtime',       0, ...
    'Comment',       'stat_2', ...
    'test_type',     'ttest_onesample', ...  % One-sample Student's t-test    X~N(m,s)t = mean(X) ./ std(X) .* sqrt(n)      df=n-1
    'tail',          'one+');  % One-tailed (+)

%% BS, Paired t-test
idx = find(contains(sFiles_A, 'ec1127')==1);

sFiles_A = sFiles_name(idx_3_smooth); sFiles_A(idx) = [];
sFiles_B = sFiles_name(idx_2_smooth);

% Process: Perm t-test paired [-0.300s,2.000s]          H0:(A=B), H1:(A<>B)
sFiles = bst_process('CallProcess', 'process_test_permutation2p', sFiles_A, sFiles_B, ...
    'timewindow',     [-0.3, 2], ...
    'scoutsel',       {}, ...
    'scoutfunc',      1, ...  % Mean
    'isnorm',         0, ...
    'avgtime',        0, ...
    'iszerobad',      1, ...
    'Comment',        '', ...
    'test_type',      'ttest_paired', ...  % Paired Student's t-test T = mean(A-B) / std(A-B) * sqrt(n)
    'randomizations', 1000, ...
    'tail',           'two');  % Two-tailed

% Process: Perm t-test paired [-1000ms,2000ms]          H0:(A=B), H1:(A<>B)
% sFiles = bst_process('CallProcess', 'process_test_permutation2p', sFiles_A, sFiles_B, ...
%     'scoutsel',       {}, ...
%     'scoutfunc',      1, ...  % Mean
%     'isnorm',         0, ...
%     'avgtime',        0, ...
%     'iszerobad',      1, ...
%     'test_type',      'ttest_paired', ...  % Paired Student's t-test T = mean(A-B) / std(A-B) * sqrt(n)
%     'randomizations', 100, ...
%     'tail',           'two');  % Two-tailed

%     'timewindow',     toi, ...

%     'Comment',        ['Perm t-test paired', ' [', num2str(toi(1)),',', num2str(toi(2)),']ms: ', stag1, '-', stag2], ...

%% FT, Paired t-test
bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles_A, sFiles_B, ...
    'timewindow',     [1, 1], ...
    'scoutsel',       {}, ...
    'scoutfunc',      1, ...  % Mean
    'isabs',          0, ...
    'avgtime',        0, ...
    'randomizations', 1000, ...
    'statistictype',  2, ...  % Paired t-test
    'tail',           'two', ...  % One-tailed (+)
    'correctiontype', 5, ...  % max
    'minnbchan',      0, ...
    'clusteralpha',   0.05);

%     'timewindow', toi, ...
%     'Comment', ['FT ttest, ', ' [', num2str(toi(1)),',', num2str(toi(2)),']ms: ', stag1, '-', stag2], ...

%% Source modelling, DICS-BF
% addpath('/data/MEG/Vahab/Github/MCW_MEGlab/Projects/ECP/SM/Surface_based/progression_analysis/function');

% disp('1: 1st');
% disp('2: 2nd');
% disp('3: 3rd');
% disp('4: 4th');
% disp('5: 5th');
% anal_sel = input(':');
%
% func = 'process_ft_sourceanalysis_DICS_BF_5intervals';
%
% datalog = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/datalog';
%
% switch anal_sel
%     case 1
%         datatag = '1st'; t_sel = 1;
%     case 2
%         datatag = '2nd'; t_sel = 2;
%     case 3
%         datatag = '3rd'; t_sel = 3;
%     case 4
%         datatag = '4th'; t_sel = 4;
%     case 5
%         datatag = '5th'; t_sel = 5;
% end
%
%
% for ii = 1:length(BS_chan.datafile)
%     [a, ~] = fileparts(BS_chan.datafile{ii});[c,d] = fileparts(a); [e,f] = fileparts(c);
%     tkz = tokenize(d,'_');
%     ee = [tkz{1},'_',tkz{2},'_',tkz{3}];
%     cd(a);
%     dd1 = rdir('./results_dics*.mat');
%
%     %     if any(strcmp(HCs,tkz{2}))
%     disp('======')
%     disp(ee)
%     %         pause,
%     run = [];
%     if ~isempty(dd1)
%         for jj=1:length(dd1)
%             tmp = load(dd1(jj).name);
%             disp(tmp.Comment);
%             run(jj) = contains(tmp.Comment, datatag);
%         end
%     end
%     runval  = isempty(find(run==1, 1));
%
%     %         if runval && ~exist(fullfile(datalog, [ee,'_', num2str(t_sel),'.mat']),'file')
%     disp(ii)
%     db_reload_studies(ii, 1);
%     %             dd = rdir('./*_IC*.mat');
%     dd = rdir('./data*.mat');
%
%     sFiles1 = [];
%     clear d_tag
%     for jj=1:length(dd)
%         sFiles1{jj} = fullfile(f,d,dd(jj).name);
%         tkz = tokenize(dd(jj).name, '_');
%         d_tag(jj) = str2num(tkz{2});
%     end
%
%     sFiles_2 = sFiles1(d_tag == 2); sFiles_3 = sFiles1(d_tag == 3);
%
%
%     w1 = 0.4; l = 1; ov = l.*1; j=1; wi=[];
%     %         w1 = 0.4; l = 0.8; ov = l.*0.2; j=1; wi=[];
%     while w1+l < 2
%         wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
%     end
%
%     for jj = 1:length(wi)
%
%         sFiles = bst_process('CallProcess', 'process_ft_sourceanalysis_dics', sFiles_2, [], ...
%             'sensortype', 'MEG', ...  % MEG
%             'poststim',   wi(jj,:), ...
%             'baseline',   [-0.3, -0.0005], ...
%             'foi',        20, ...
%             'tpr',        4, ...
%             'method',     'subtraction', ...  % Subtraction (post-pre)
%             'erds',       'erd', ...  % ERD
%             'effect',     'abs', ...  % abs
%             'Comment', [ee,'_2_', num2str(t_sel),'.mat'], ...
%             'maxfreq',    40, ...
%             'showtfr',    1);
%
%
%         sFiles = bst_process('CallProcess', 'process_ft_sourceanalysis_dics', sFiles_3, [], ...
%             'sensortype', 'MEG', ...  % MEG
%             'poststim',   wi(jj,:), ...
%             'baseline',   [-0.3, -0.0005], ...
%             'foi',        20, ...
%             'tpr',        4, ...
%             'method',     'subtraction', ...  % Subtraction (post-pre)
%             'erds',       'erd', ...  % ERD
%             'effect',     'abs', ...  % abs
%             'Comment', [ee,'_3_', num2str(t_sel),'.mat'], ...
%             'maxfreq',    40, ...
%             'showtfr',    1);
%
%     end
%
%     %                         pause,
%     save(fullfile(datalog, [ee,'_Sybl_', num2str(t_sel),'.mat']),'t_sel'); % prevent conflict with other running scripts
%     bst_process('CallProcess', func, sFiles_2, [], ...
%         'method',     'dics', ...  % DICS beamformer
%         'sensortype', 'MEG', ...
%         't_sel', t_sel, 'progressbar', 0);  % MEG
%     disp('done');
%     delete(fullfile(datalog, [ee,'_', num2str(t_sel),'.mat'])); % delete tmp file
%     clc,
%     %         end
%     %     end
% end
%
% subj_all1 = unique(subj_all);
