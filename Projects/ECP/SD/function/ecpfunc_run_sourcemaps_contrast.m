function ecpfunc_run_sourcemaps_contrast(cfg_main)

protocol = cfg_main.protocol;
% data_info_dir = cfg_main.datadir;
BS_data_dir = cfg_main.BS_data_dir;

%%
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

sub_all1 = setdiff(unq_bs_subj,no_anat);

%%
cd(BS_data_dir)
d = dir('./Group_analysis/*SD*');

for j=1:length(d)
    cd(BS_data_dir)
    
    cd(fullfile('Group_analysis', d(j).name))
    
    dd = rdir(fullfile('./results*ssmooth.mat'));
    dd1 = rdir(fullfile('./results_abs*.mat'));
    
    if length(dd) == 2 && isempty(dd1)
        
        Comment = [];
        for jj=1:length(dd)
            disp([num2str(jj), '/' , num2str(length(dd))])

            tmp  = load(dd(jj).name);
            Comment{jj} = tmp.Comment;
        end
        disp(Comment)
        if find(contains(Comment, 'Avg: 3')==1) && find(contains(Comment, 'Avg: 2')==1)
            
            idx_3 = contains(Comment, 'Avg: 3')==1;
            idx_2 = contains(Comment, 'Avg: 2')==1;
            
            [pathstr, name, ~] = fileparts(dd(jj).folder);
            [~, name2, ~] = fileparts(pathstr);
            
            sFiles_3 = fullfile(name2, name, dd(idx_3).name);
            sFiles_2 = fullfile(name2, name, dd(idx_2).name);
            
            % Process: Difference: A-B, abs
            bst_process('CallProcess', 'process_diff_ab', sFiles_3, sFiles_2, ...
                'source_abs', 1);      
        end
    end
end

%%
cd(BS_data_dir)
dd = rdir(('./Group_analysis/*SD*/results_abs*.mat'));

clc
Comment = [];
for j=1:length(dd)
    disp([num2str(j), '/' , num2str(length(dd))])   
    tmp  = load(dd(j).name);
    Comment{j} = tmp.Comment;
    sub{j} = tmp.Comment(9:14);
    disp(Comment{j})
    disp(sub{j})
end

%%
unq_sub = unique(sub);

cd(BS_data_dir)

for i=1:length(unq_sub)
    idx = find(contains(sub,unq_sub{i})==1);
    if length(idx) == 2
        Comment1 = [];
        for j=1:length(idx)
            tmp = load(dd(idx(j)).name);
            Comment1{j} = tmp.Comment;
        end
        disp(Comment1')
        
        sFiles_run1 = dd(idx(1)).name(3:end);
        sFiles_run2 = dd(idx(2)).name(3:end);
        
        disp(sFiles_run1)
        disp(sFiles_run2)
        
        sFiles = {sFiles_run1; sFiles_run2};
        
        % Process: Average: Everything
        bst_process('CallProcess', 'process_average', sFiles, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...
            'scalenormalized', 0);
    end
end