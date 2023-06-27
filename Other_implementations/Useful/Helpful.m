
%% create newfolder
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

%% SCP copy
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

%% Asking input
disp('extract cooditnes? yes=1, no=2'); ask_extractcoor = input('');
if ask_extractcoor ==1
end

%% Intersect (2 strings)
[C,ia,ib] = intersect(sort_sub_name, sort_sub_name2, 'stable');

%% Lowercase and uppercase 
lower(sub_name2);

%% alphabetically sorting cells
[sort_sub_name2,~] = sort(lower(sub_name2));

%% check if a variable exist in a cell matlab
isfield(a, 'History')

%% Check if a folder exists
isfolder(['@rawec',sub_sel, '_SD_', run_sel, '_raw'])

%% FT Progress 
clc
ft_progress('init', 'text',     'please wait ...');
for i=1:42
    ft_progress(i/42, 'Processing event %d from %d', i, 42);
    pause(0.1);
end
ft_progress('close');

%% Normalize values
tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:)));

%% Java comments
javaclasspath('-remove','/usr/local/MATLAB_Tools/brainstorm3/java/brainstorm.jar')
javarmpath('/usr/local/MATLAB_Tools/brainstorm3/java/brainstorm.jar');
pp = javaclasspath('-dynamic');
for i=1:length(pp)
    javarmpath(pp{i})
end
javaclasspath('-dynamic')
javaclasspath

%% Remove tool from the path
rmpath(spm_path);

%% Table
t1 = table(m_LI_sub'); t1.Properties.VariableNames{'Var1'} = 'LI';

%% Create a folde

folderPath = 'path/to/folder';

if ~exist(folderPath, 'dir')
    mkdir(folderPath);
    disp('Folder created successfully.');
else
    disp('Folder already exists.');
end

