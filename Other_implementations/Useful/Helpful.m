
%% Create newfolder
if exist(savepath, 'file') == 0, mkdir(savepath), end

%% Counting folders
all_files = dir;
all_dir = all_files([all_files(:).isdir]);
num_dir = numel(all_dir);

if numel(dir(path)) <=2
    % folder is empty
else
    % folder isn't empty
end

%% SCP copy (between/within servers)
% - copy data between servers
command = 'scp -r /MEG_data/MRI_database/epilepsy/DEJESUS_Elias_T1_AxBravo/nii vyoussofzadeh@squiggles.rcc.mcw.edu:/data/MEG/Vahab/';
% command = 'scp -r /data/MEG/Vahab/Data_clinical/CAT_analysis/DAMICO_Madalynn_AxT1 vyoussof@192.227.57.19:/MEG_data/MRI_database/epilepsy/DAMICO_Madalynn_AxT1/CAT_10_07_2022';
system(command)

% - copy data between folders
scp -r /MEG_data/Software/CAT12/cat12 /usr/local/MATLAB_Tools/spm12/toolbox

%% Regression
tbl = table(X(:,idx(1)),y1');

mdl = fitlm(tbl,'linear');
plot(mdl)

%% Asking for input
disp('extract cooditnes? yes=1, no=2'); ask_extractcoor = input('');
if ask_extractcoor ==1
end

%% Intersect (2 strings)
[C,ia,ib] = intersect(sort_sub_name, sort_sub_name2, 'stable');

%% Lowercase and uppercase of a string cell
lower(sub_name2);

%% alphabetically sorting cells
[sort_sub_name2,~]=sort(lower(sub_name2));

%% File exists in a folder
if ~exist(fullfile(savepath, tkz{end}),'file')
end

%% Copyfile 
copyfile(d_tar.name,savedir),

%% remove space from string matlab