function S_data = ecpfunc_read_sourcemaps(cfg)

protocol = cfg.protocol;

%%

% db_reload_database('current',1)
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
    'EC1092'
    'EC1094'
    'EC1096'
    'EC1110'
    'EC1111'
    'EC1112'
    'EC1141'
    'EC1153'
    'EC1162'
    'EC1090'};

sub_all1 = setdiff(unq_bs_subj,no_anat);

%%
data_info_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas';

cd(data_info_dir)
if exist(fullfile(data_info_dir,'comments_subject.mat'),'file') == 2
    load(fullfile(data_info_dir,'comments_subject.mat')),
else
    subj_del = [];
    
    cd(BS_data_dir)
    dd = rdir(fullfile('./Group_analysis/1_LCMV_Subjects/results_average*.mat'));
    for jj=1:length(dd), disp([num2str(jj),':',dd(jj).name]); end
    
    sFiles_name = [];
    for jj=1:length(dd)
        sFiles_name{jj} = fullfile(dd(jj).name(3:end));
    end
    Comment = []; k=1; kk=1;
    need_correction = [];
    for jj=1:length(sFiles_name)
        disp([num2str(jj), '/' , num2str(length(sFiles_name))])
        cd(BS_data_dir)
        tmp  = load(sFiles_name{jj});
        Comment{jj} = tmp.Comment;
        if length(tmp.Time) < 4000 && ~contains(Comment{jj}, 'Avg')
            disp(length(tmp.Time))
            need_correction(k) = jj;
            k=k+1;
        end
        if contains(Comment{jj}, 'Avg:')
            Comment_sel{kk} = sFiles_name{jj};
            kk=kk+1;
            disp(Comment{jj})
        end
    end
    
    removing_sFiles  = sFiles_name(need_correction);
    for i=1:length(removing_sFiles)
        %         removing_sFiles{i};
        [path, ~] = fileparts(removing_sFiles{i});
        if exist(path,'dir') == 7
            rmdir(fullfile(BS_data_dir, path),'s')
        end
    end
    
    idx_anim = find(contains(Comment, 'Anim_')==1);
    idx_symb = find(contains(Comment, 'Symbol_')==1);
    
    sFiles_3 = sFiles_name(idx_anim);
    sFiles_2 = sFiles_name(idx_symb);
    
    save(fullfile(data_info_dir,'comments_subject.mat'),'Comment','sFiles_name','sFiles_3','sFiles_2'),
end

%%
cd(data_info_dir)
if exist(fullfile(data_info_dir,'subjects_ID.mat'),'file') == 2
    load(fullfile(data_info_dir,'subjects_ID.mat')),
else
    cd(BS_data_dir)
    L = length(sFiles_3); k = 1;
    clear subjs_3
    for i=1:length(sFiles_3)
        %     disp([num2str(i), '/' , num2str(length(sFiles_3))])
        tmp = load(sFiles_3{i});
        idx = strfind(tmp.Comment,'_');
        subjs_3{k} = tmp.Comment(idx+1:idx+6);
        disp([subjs_3{k}, ': ', num2str(i), '/' , num2str(length(sFiles_3))])
        k=1+k;
    end
    unq_bs_subj_3 = unique(subjs_3);
    
    L = length(sFiles_2); k = 1;
    clear subjs_2
    for i=1:length(sFiles_2)
        %     disp([num2str(i), '/' , num2str(length(sFiles_3))])
        tmp = load(sFiles_2{i});
        idx = strfind(tmp.Comment,'_');
        subjs_2{k} = tmp.Comment(idx+1:idx+6);
        disp([subjs_2{k}, ': ', num2str(i), '/' , num2str(length(sFiles_2))])
        k=1+k;
    end
    unq_bs_subj_2 = unique(subjs_2);
    
    save(fullfile(data_info_dir,'subjects_ID.mat'),'subjs_2','subjs_3'),
end

%%
S_data = [];
S_data.sFiles_3 = sFiles_3;
S_data.sFiles_2 = sFiles_2;
S_data.subjs_3 = subjs_3;
S_data.subjs_2 = subjs_2;

