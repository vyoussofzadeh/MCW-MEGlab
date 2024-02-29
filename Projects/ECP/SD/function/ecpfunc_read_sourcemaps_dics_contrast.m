function S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg_main)

protocol = cfg_main.protocol;
data_info_dir = cfg_main.datadir;
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

% sub_all1 = setdiff(unq_bs_subj,no_anat);

%%
cd(data_info_dir)

subj_del = [];

cd(BS_data_dir)
dd = rdir(['./Group_analysis/', cfg_main.datatag, '/results*.mat']);
d_sub = rdir(['./Group_analysis/', cfg_main.datatag, '/results*_mean.mat']);


% Extract names
dd_names = {dd.name};
d_sub_names = {d_sub.name};

% Find and remove matching entries
[~, idx_to_remove] = ismember(dd_names, d_sub_names);
dd(idx_to_remove > 0) = [];


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

idx_anim_symb = contains(Comment, 'wDICS_')==1;
sFiles_32 = sFiles_name(idx_anim_symb);


%%
cd(data_info_dir)


cd(BS_data_dir)
L = length(sFiles_32); k = 1;
clear subjs
for i=1:length(sFiles_32)
    tmp = load(sFiles_32{i});
    idx = strfind(tmp.Comment,'_');
    subjs{k} = tmp.Comment(idx(1)-6:idx(1)-1);
    disp([subjs{k}, ': ', num2str(i), '/' , num2str(length(sFiles_32))])
    k=1+k;
end
unq_bs_subj = unique(subjs);


%%
S_data = [];
S_data.sFiles_32 = sFiles_32;
S_data.subjs = subjs;

