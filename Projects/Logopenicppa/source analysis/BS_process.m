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
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Logopenicppa/functions');

%%
flags = [];
flags.ica = 0;

%%
% adding BS path
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3'; %BS-2021
addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

BS_dir = '/data/MEG/Research/logopenicppa/BS_process/Logopenicppa_CRM_zscore';
BS_data_dir = fullfile(BS_dir,'data');
protocol = fullfile(BS_dir, 'data/protocol.mat');

%% Loading up raw data
tag = 'CRM'; %- Semantic decision

cfg = [];
cfg.indir = indir;
cfg.tag = tag;
datafile = do_datalookup(cfg);

%%
db_reload_database('current',1)
load(protocol);
Subj_bs = ProtocolSubjects.Subject;

datafile = [];
subjs_bs = [];

L = length(Subj_bs);
k = 1;
clear subjs_bs
for i=1:L
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
if flags.ica == 1
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
        
        dd = rdir(fullfile(pathstr,'data_3*_trl.mat'));
        if isempty(dd)
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

%% Source, lcmv
d = rdir(fullfile(BS_data_dir,'/*/*CRM*clean_low/channel_vectorview306_acc1.mat'));

ft_progress('init', 'text', 'please wait ...');
for i=1:length(d)
    
    ft_progress(i/length(d), 'Processing event %d from %d', i, length(d));
    
    pause(0.1);
    [pathstr, name] = fileparts(d(i).name);
    cd(pathstr)
    dd_results = rdir(fullfile(pathstr,'results*PNAI*.mat'));
    results = [];
    for j=1:length(dd_results)
        cmt = load(dd_results(j).name);
        results{j} = cmt.Comment;
    end
    
    dd = rdir(fullfile(pathstr,'data_1*_trial*.mat'));
    if isempty(dd)
        dd = rdir(fullfile(pathstr,'data_1*_trial*.mat'));
        if isempty(dd)
            dd = rdir(fullfile(pathstr,'data_1*_trial*.mat'));
        end
    end
    
    if isempty(results)
        flag = 1;
    else
        flag = 0;
    end
    
    if ~isempty(dd) && flag == 1
        pause,
        sFiles = [];
        for j=1:length(dd)
            [pathstr, name] = fileparts(dd(j).name);
            [pathstr2, name2] = fileparts(pathstr);
            [~, name3] = fileparts(pathstr2);
            sFiles{j} = fullfile(name3, name2, [name, '.mat']);
        end
        cfg = []; cfg.comment = 'CRM'; bs_lcmv(cfg,sFiles)
    end
end
ft_progress('close');


%% Z-score normalization (it could be done in combination with bs_lcmv)

d = rdir(fullfile(BS_data_dir,'/*/*CRM*clean_low/channel_vectorview306_acc1.mat'));


ft_progress('init', 'text', 'please wait ...');
for i=1:length(d)
    
    ft_progress(i/length(d), 'Processing event %d from %d', i, length(d));
    pause(0.1);
    
    [pathstr, name] = fileparts(d(i).name);
    cd(pathstr)
    dd_results = rdir(fullfile(pathstr,'results*PNAI*KERNEL*.mat'));
    results = [];
    for j=1:length(dd_results)
        cmt = load(dd_results(j).name);
        results{j} = cmt.Comment;
    end
    
    if ~isempty(results)
        flag = 1;
    else
        flag = 0;
    end
    
    dd = rdir(fullfile(pathstr,'data_1_average*.mat'));
    ddd = rdir(fullfile(pathstr,'*_zscore*'));
    
    if flag == 1 && ~isempty(dd) && isempty(ddd)
                
        [pathstr, name] = fileparts(dd_results.name); [pathstr2, name2] = fileparts(pathstr); [~, name3] = fileparts(pathstr2);
        sFiles_source = fullfile(name3, name2, [name, '.mat']);
        
        [pathstr_avg, name_avg] = fileparts(dd.name); [pathstr2_avg1, name_avg2] = fileparts(pathstr_avg); [~, name_avg3] = fileparts(pathstr2_avg1);
        sFiles_avg = fullfile(name_avg3, name_avg2, [name_avg, '.mat']);
        
        sFiles = {['link|',sFiles_source, '|',sFiles_avg]};
        
        disp(sFiles)
        bst_process('CallProcess', 'process_baseline_norm', sFiles, [], ...
            'baseline',   [-0.3, -0.0005], ...
            'source_abs', 0, ...
            'method',     'zscore', ...  % Z-score transformation:    x_std = (x - &mu;) / &sigma;
            'overwrite',  1);       
    end
end
ft_progress('close');

%% Avg and Project to default template
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

for ii = 1:length(subj)
    
    cd(fullfile(BS_data_dir,subj{ii}))
    dd = rdir(['./*',tag,'*/results_PNAI*zscore.mat']);
    for jj=1:length(dd), disp([num2str(jj),':',dd(jj).name]); end
    
    sFiles1 = []; sFiles_name = [];
    for jj=1:length(dd)
        sFiles1{jj} = fullfile(subj{ii},dd(jj).name);
        [aa,bb] = fileparts(dd(jj).name);
        sFiles_name{jj} = aa;
    end
    
    [pathstr, ~, ~] = fileparts(sFiles1{1}); 
    [pathstr, name, ext] = fileparts(pathstr);
    
    ddd = rdir('./@intra/results*.mat');
    
    for j = 1:length(ddd)
        disp(tmp.Comment)
        tmp = load(ddd(j).name);
        if contains(tmp.Comment,['Avg: ',pathstr,'_zscored'])            
            flag = 0; break,
        else
            flag = 1;
        end
    end
    
    if length(sFiles_name) == 2 && flag == 1
        
        pause,
        % Process: Average: Everything
         bst_process('CallProcess', 'process_average', sFiles1, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...
            'Comment', [pathstr, '_zscored'], ...
            'scalenormalized', 0);
    end
end

%% Default template
for ii = 1:length(subj)
    
    cd(fullfile(BS_data_dir,subj{ii}))
    ddd = rdir('./@intra/results*.mat');

    for j = 1:length(ddd)
        tmp = load(ddd(j).name);
        disp(tmp.Comment)
        if contains(tmp.Comment,'_zscored')            
            flag = 1; k = j;break,
        else
            flag = 0;
        end
    end
    if length(sFiles_name) == 2 && flag == 1
        sFiles1 = fullfile(subj(ii), ddd(k).name);
        bst_project_sources(sFiles1, destSurfFile, 0, 1);
    end
end

