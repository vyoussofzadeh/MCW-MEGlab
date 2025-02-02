function S_data = ecpfunc_read_sourcemaps2(cfg_main)

data_info_dir = cfg_main.datadir;
BS_data_dir = cfg_main.BS_data_dir;


%%
cd(BS_data_dir)
% dd = rdir(fullfile('./Group_analysis/ec*/results*abs_ssmooth.mat'));

dd = rdir(cfg_main.datamask);
for jj=1:length(dd), disp([num2str(jj),':',dd(jj).name]); end

sFiles_name = [];
for jj=1:length(dd)
    sFiles_name{jj} = fullfile(dd(jj).name(3:end));
end

Comment = []; kk=1;
for jj=1:length(sFiles_name)
    disp([num2str(jj), '/' , num2str(length(sFiles_name))]);
    cd(BS_data_dir)
    tmp  = load(sFiles_name{jj});
    
    % Extract run number from filename
%     run_str = strfind(sFiles_name{jj}, '_run'); sFiles_name{jj}(run_str+4);
    Comment{jj} = tmp.Comment;
%     Comment{jj} = [tmp.Comment, ' Run: ', sFiles_name{jj}(run_str+4)]; % Append run_number to Comment
    if contains(Comment{jj}, 'Avg:')
        Comment_sel{kk} = sFiles_name{jj};
        kk=kk+1;
        disp(Comment{jj})
    end
end

idx_anim = contains(Comment, '3 (')==1;
idx_symb = contains(Comment, '2 (')==1;

sFiles_3 = sFiles_name(idx_anim);
sFiles_2 = sFiles_name(idx_symb);

%%
cd(data_info_dir)
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

%%
S_data = [];
S_data.sFiles_3 = sFiles_3;
S_data.sFiles_2 = sFiles_2;
S_data.subjs_3 = subjs_3;
S_data.subjs_2 = subjs_2;


