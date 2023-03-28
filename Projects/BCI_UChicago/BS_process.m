%% BCI data, Medical College of Wisconsin - Uni Chicago

% Script: BS Process (preprocessing, source analysis)
% Project: BCI
% Writtern by: Vahab YoussofZadeh
% Update: 03/26/2023

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
cd('/MEG_data/Research_studies/BCI_Uni_of_Chicago')

indir = '/MEG_data/Research_studies/BCI_Uni_of_Chicago/bci03_bci03';
outdir = '/MEG_data/Research_studies/BCI_Uni_of_Chicago/Results';

% ft_path = '/usr/local/MATLAB_Tools/fieldtrip_20190419';
ft_path = '/MEG_data/Software/Fieldtrip/latest/fieldtrip-master';
addpath(ft_path); ft_defaults

%%
ft_func = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/FT_fucntions';
addpath(genpath(fullfile(ft_func, '/External')));
addpath(fullfile(ft_func, '/functions_new'));

addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/BCI_UChicago/functions')


%% Lookup raw data
cfg = [];
cfg.indir = indir;
cfg.tag = tag;
datafile = do_datalookup(cfg);

%%
% adding BS path
bs_path = '/MEG_data/Software/Brainstorm/Brainstorm_2022/brainstorm3';

% bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3'; %BS-2021
addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

BS_dir = '/MEG_data/Research_studies/BCI_Uni_of_Chicago/BS_database';
BS_data_dir = fullfile(BS_dir,'data');
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
disp(unq_bs_subj)

sub_all1 = unq_bs_subj;

%% Est. head model
cd(BS_data_dir)
d = rdir(fullfile(BS_data_dir,'/*/Run*_band/channel_vectorview306_acc1.mat'));

OPTIONS = [];
OPTIONS.Comment = 'Overlapping spheres';
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

%% Clean-up (ICA)
addpath(ft_path);
ft_defaults

d = rdir(fullfile(BS_data_dir,'/*/Run*_band/channel_vectorview306_acc1.mat'));

for i=1:length(d)
    
    [pathstr, name] = fileparts(d(i).name);
    cd(pathstr)
    
    dd_results = rdir(fullfile(pathstr,'results*.mat'));
    results_comments = [];
    for jj=1:length(dd_results)
        tmp = load(dd_results(jj).name);
        results_comments{jj} = tmp.Comment(1);
    end
    
    dd = rdir(fullfile(pathstr,'data_1*_ic.mat'));
    if isempty(dd)
        dd_1 = rdir(fullfile(pathstr,'data_1*.mat'));
        sFiles = [];
        for j=1:length(dd_1)
            [pathstr, name] = fileparts(dd_1(j).name);
            [pathstr2, name2] = fileparts(pathstr);
            [pathstr3, name3] = fileparts(pathstr2);
            sFiles{j} = fullfile(name3, name2, [name, '.mat']);
        end
        disp(sFiles(i))
        idx = strfind(name,'_');
        if ~isempty(results_comments)
            runidx = contains(results_comments, name(idx(1)+1));
        else
            runidx = [];
        end
        if isempty(runidx)
            %             aks_ica = input('run ICA (yes=1)?');
            aks_ica = 1;
            if aks_ica ==1
                bs_preprocess_ICA_trials(sFiles),
            else
                disp('done')
            end
        end
    end
end
db_reload_database('current',1)

%% Clean-up (reject trials)
d = rdir(fullfile(BS_data_dir,'/*/Run*_band/channel_vectorview306_acc1.mat'));

for i=1:length(d)
    
    [pathstr, name] = fileparts(d(i).name);
    cd(pathstr)
    
    dd_results = rdir(fullfile(pathstr,'results*.mat'));
    results_comments = [];
    for jj=1:length(dd_results)
        tmp = load(dd_results(jj).name);
        results_comments{jj} = tmp.Comment(1);
    end
    
    dd = rdir(fullfile(pathstr,'data_1*_trl.mat'));
    if isempty(dd)
        dd_1 = rdir(fullfile(pathstr,'data_1*_ic.mat'));
        sFiles = [];
        for j=1:length(dd_1)
            [pathstr, name] = fileparts(dd_1(j).name);
            [pathstr2, name2] = fileparts(pathstr);
            [pathstr3, name3] = fileparts(pathstr2);
            sFiles{j} = fullfile(name3, name2, [name, '.mat']);
        end
        disp(sFiles(i))
        idx = strfind(name,'_');
        if ~isempty(results_comments)
            runidx = contains(results_comments, name(idx(1)+1));
        else
            runidx = [];
        end
        if isempty(runidx)
            %             aks_ica = input('run rej trls (yes=1)?');
            aks_ica = 1;
            if aks_ica ==1
                bs_preprocess_reject_trials(sFiles),
            else
                disp('done')
            end
        end
    end
end
db_reload_database('current',1)

%% Source modelling, DICS-BF
d = rdir(fullfile(BS_data_dir,'/*/Run*_band/channel_vectorview306_acc1.mat'));


for i=1:length(d)
    
    [pathstr, name] = fileparts(d(i).name);
    
    dd_results = rdir(fullfile(pathstr,'*subtraction_MEG*.mat'));
    results = [];
    for j=1:length(dd_results)
        cmt = load(dd_results(j).name);
        results{j} = cmt.Comment;
    end
    
    dd = rdir(fullfile(pathstr,'data_1*ic_trl.mat'));
    if isempty(dd)
        dd = rdir(fullfile(pathstr,'data_1_ic*.mat'));
        if isempty(dd)
            dd = rdir(fullfile(pathstr,'data_1_trial*.mat'));
        end
    end
    
    if isempty(results)
        flag = 1;
    elseif isempty(find(contains(results, 'subtraction')==1, 1))
        flag = 1;
    else
        flag = 0;
    end
    
    if ~isempty(dd) && flag == 1
        sFiles = [];
        for j=1:length(dd)
            [pathstr, name] = fileparts(dd(j).name);
            [pathstr2, name2] = fileparts(pathstr);
            [~, name3] = fileparts(pathstr2);
            sFiles{j} = fullfile(name3, name2, [name, '.mat']);
        end
        
        sFiles = bst_process('CallProcess', 'process_ft_sourceanalysis_dics', sFiles, [], ...
            'sensortype', 'MEG', ...  % MEG
            'poststim',   [4.274358645e-15, 1], ...
            'baseline',   [-0.3, -0.001], ...
            'foi',        20, ...
            'tpr',        4, ...
            'method',     'subtraction', ...  % Subtraction (post-pre)
            'erds',       'erd', ...  % ERD
            'effect',     'abs', ...  % abs
            'maxfreq',    40, ...
            'showtfr',    1);
    end
end
close all

%% Export maps - individual tasks
d = rdir(fullfile(BS_data_dir,'/*/Run*_band/channel_vectorview306_acc1.mat'));
for i=1:length(d)
    
    [pathstr, name] = fileparts(d(i).name);
    
    dd_results = rdir(fullfile(pathstr,'*subtraction_MEG*.mat'));
    results = [];
    for j=1:length(dd_results)
        cmt = load(dd_results(j).name);
        results{j} = cmt.Comment;
    end
    
    FileName  = dd_results.name;
    [sStudy, iStudy, iItem] = bst_get('AnyFile', FileName);
    ProtocolInfo = bst_get('ProtocolInfo');
    [FileName, FileType, isAnatomy] = file_fullpath(FileName);
    filePath = ProtocolInfo.STUDIES;
    fileBase = file_win2unix(strrep(FileName, filePath, ''));
    tkz = tokenize(fileBase,'/'); tkz1 = tokenize(tkz{2},'_');
    
    cfg = [];
    cfg.svdir = outdir;
    cfg.BSpath = filePath;
    cfg.fname = fileBase;
    cfg.side_sel = 2;
    cfg.svname = [tkz1{1}, '_', tkz1{2}];
    cfg.title = [tkz1{1}, '_', tkz1{2}];
    do_export_images_BS(cfg)
    cd(outdir)
end


%% Export maps - avg tasks
d = rdir(fullfile(BS_data_dir,'/*/@intra/results*.mat'));

for i=1:length(d)
    
    [pathstr, name] = fileparts(d(i).name);
    
    FileName  = d(i).name;
    [sStudy, iStudy, iItem] = bst_get('AnyFile', FileName);
    ProtocolInfo = bst_get('ProtocolInfo');
    [FileName, FileType, isAnatomy] = file_fullpath(FileName);
    filePath = ProtocolInfo.STUDIES;
    fileBase = file_win2unix(strrep(FileName, filePath, ''));
    
    tmp = load(FileName);
    tkz = tokenize(tmp.Comment,'_');
    
    disp(tkz{end})
    
    cfg = [];
    cfg.svdir = outdir;
    cfg.BSpath = filePath;
    cfg.fname = fileBase;
    cfg.side_sel = 2;
    cfg.svname = tkz{end};
    do_export_images_BS(cfg)
    cd(outdir)
end
