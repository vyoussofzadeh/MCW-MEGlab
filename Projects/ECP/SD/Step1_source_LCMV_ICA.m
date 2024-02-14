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

%% Source, lcmv
% d1 = rdir(fullfile(BS_data_dir,'/*/ec*raw_low_clean/channel_vectorview306_acc1.mat'));
% d2 = rdir(fullfile(BS_data_dir,'/*/ec*raw_tsss_low_clean/channel_vectorview306_acc1.mat'));
% d3 = rdir(fullfile(BS_data_dir,'/*/EC*raw_low_clean/channel_vectorview306_acc1.mat'));
% d = [d1;d2;d3];

d1 = rdir(fullfile(BS_data_dir,'/*/ec*clean_low_clean/channel_vectorview306_acc1.mat'));
d2 = rdir(fullfile(BS_data_dir,'/*/ec*raw_tsss_low_clean/channel_vectorview306_acc1.mat'));
d3 = rdir(fullfile(BS_data_dir,'/*/EC*raw_low_clean/channel_vectorview306_acc1.mat'));
d4 = rdir(fullfile(BS_data_dir,'/*/ec*raw_clean_low/channel_vectorview306_acc1.mat'));
d5 = rdir(fullfile(BS_data_dir,'/*/ec*elecfix*raw_clean*/channel_vectorview306_acc1.mat'));

d = [d1;d2;d3; d4; d5];

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
%             pause
%             sFiles = bst_process('CallProcess', 'process_project_sources', sFiles2, [], ...
%                 'headmodeltype', 'surface');  % Cortex surface
            sFiles = bst_project_sources({sFiles2}, destSurfFile, 0, 1);
        end
    end
end

%% delete extra smoothed files - Optional
cd(BS_data_dir)
clc
dd = rdir(fullfile('./Group_analysis/*_clean/results_*.mat'));

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

%% delete extra avg file - Optional
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

%% delete short epochs - Optional
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
dd = rdir(fullfile('./Group_analysis/ec*clean/results_*.mat'));
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
% ft_progress('close');

%% Contrast data conditions 3 - 2
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/function')
datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';

cfg = [];
cfg.protocol = protocol;
cfg.datadir = datadir;s
cfg.BS_data_dir = BS_data_dir;
S_data = ecpfunc_read_sourcemaps(cfg);

bst_process('CallProcess', 'process_diff_ab', S_data.sFiles_3, S_data1.sFiles_2, ...
    'source_abs', 1);

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


ft_progress('init', 'text',     'please wait ...');
sFiles_name = [];
for jj=1:length(dd)
    ft_progress(jj/length(dd), 'Processing source files %d from %d', jj, length(dd));
    sFiles_name{jj} = fullfile(dd(jj).name(3:end));
    pause(0.1);
end
ft_progress('close');

ft_progress('init', 'text',     'please wait ...');
Comment = [];
for jj=1:length(sFiles_name)
    ft_progress(jj/length(dd), 'Processing source files %d from %d', jj, length(dd));
    cd(BS_data_dir)
    tmp  = load(sFiles_name{jj});
    Comment{jj} = tmp.Comment;
    pause(0.1);
end
ft_progress('close');

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

%%
