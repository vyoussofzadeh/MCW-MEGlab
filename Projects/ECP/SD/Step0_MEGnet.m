

cd('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD')
addpath('./run')
Run_setpath
addpath('./data')
addpath('./run')
% addpath(genpath('./functions'))

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/helper')

%%
% cd('/data/MEG/Research/MEGnet')
% command = 'ICA.py -filename $(pwd)/sss/ec1002_SD_run1_raw.fif -results_dir $(pwd)/outputs -line_freq  60';
% system(command)

%% tSSS files (SD task)
% /group/jbinder/ECP/MEG/MEG_Work/EC1114/tSSS/ec1114_SD_run1_raw.fif

% /data/MEG/Research/ECP/Semantic_Decision/Raw_data

addpath(genpath('/data/MEG/Research/Rupesh/Scripts'))

%%
datadir = '/group/jbinder/ECP/MEG/MEG_Work';

%% Run MEGnet analysis
rawfiflist = rdir(fullfile(datadir,'/*/tSSS/*SD*raw.fif'));


clc
for i=1:length(rawfiflist)
    rawfiflist_name = rawfiflist(i).name;
    idx = strfind(rawfiflist_name,'/');
    [pathstr, name, ext] = fileparts(rawfiflist_name);
    ecp_file = rawfiflist_name(idx(end)+1:end);
    cd(pathstr),
    cd ..
    idx = strfind(ecp_file,'_run');
    if isempty(idx)
        idx = strfind(ecp_file,'_Run');
    end
    run_num = ecp_file(idx(1)+1:idx(1)+4);
    if exist('ICA', 'file') == 0, mkdir('ICA'), end
    if isempty(rdir(['./ICA/*', run_num, '*/ica_clean.fif']))
        command = ['ICA.py -filename $(pwd)/tSSS/', ecp_file, ' -results_dir $(pwd)/ICA -line_freq  60'];
        disp(['MEGnet Processing of ', ecp_file])
        system(command)
    else
        disp(['ICA data for ', name,' exist, skipped!'])
    end
end

rawfiflist_2 = rdir(fullfile(datadir,'/*/tSSS/*SD*tsss.fif'));
clc
for i=1:length(rawfiflist_2)
    rawfiflist_name = rawfiflist_2(i).name;
    idx = strfind(rawfiflist_name,'/');
    [pathstr, name, ext] = fileparts(rawfiflist_name);
    ecp_file = rawfiflist_name(idx(end)+1:end);
    cd(pathstr),
    cd ..
    idx = strfind(ecp_file,'_run');
    if isempty(idx)
        idx = strfind(ecp_file,'_Run');
    end
    run_num = ecp_file(idx(1)+1:idx(1)+4);
    if exist('ICA', 'file') == 0, mkdir('ICA'), end
    if isempty(rdir(['./ICA/*', run_num, '*/ica_clean.fif']))
        command = ['ICA.py -filename $(pwd)/tSSS/', ecp_file, ' -results_dir $(pwd)/ICA -line_freq  60'];
        disp(['MEGnet Processing of ', ecp_file])
        system(command)
    else
        disp(['ICA data for ', name,' exist, skipped!'])
    end
end

%% Copy and rename the ICA_clean data

clc
%- listing raw data
tag = 'SD'; %- Semantic decision
d = rdir([datadir,['/*/ICA/*/ica_clean.fif']]);

clear datafolder datafile sub subj_all
datafile_fif = [];
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    sub{i} = datafile{i}(Index(7)+1:Index(8)-1);
    subj_all{i} = [num2str(i), ': ', sub{i}];
end
datafile_fif = vertcat(datafile_fif,datafile);
disp(datafile_fif')

% Rename (update naming of ICA_clean file with sub-IDs)
newnames = cell(size(datafile_fif)); % Preallocate for efficiency
for i = 1:length(datafile_fif)
    tkz = tokenize(datafile_fif{i}, '/');
    [pathstr, name, ext] = fileparts(datafile_fif{i});
    % Construct the new file name
    newnames{i} = fullfile(pathstr, [tkz{end-1}, '_', name, ext]);

    % Check if the file has already been renamed
    if ~exist(newnames{i}, 'file')
        disp(['Renaming: ', datafile_fif{i}])
%         pause,
        commandStr = ['cp ', datafile_fif{i}, ' ', newnames{i}];
        [status, cmdout] = system(commandStr, '-echo');
        if status == 0
            disp('File copied and renamed.');
        else
            disp(['Error in copying file: ', cmdout]);
        end
    else
        disp([datafile_fif{i}, ' has already been renamed!'])
    end
end

%% Copy and rename the ICA_clean data
clc
%- listing raw data
tag = 'SD'; %- Semantic decision
d = rdir([datadir,['/*/ICA/*/ica_clean.fif']]);

clear datafolder datafile sub subj_all
datafile_fif = [];
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    sub{i} = datafile{i}(Index(7)+1:Index(8)-1);
    subj_all{i} = [num2str(i), ': ', sub{i}];
end
datafile_fif = vertcat(datafile_fif,datafile);
disp(datafile_fif')

%- Rename (update naeming of ICA_clean file with sub-IDs)
clc
newnames = [];
for i=1:length(datafile_fif)
    tkz = tokenize(datafile_fif{i},'/');
    [pathstr, name] = fileparts(datafile_fif{i});
    newnames{i} = fullfile(pathstr, [tkz{end-1}(1:end), '_', tkz{end}]);
    cd(pathstr)
    if isempty(rdir(newnames{i}))
        disp(newnames{i})
        pause,
%         disp(newnames{i})
%         copyfile(datafile_fif{i},newnames{i},'f');
%         disp('copied.')
    else
%         disp(datafile_fif{i})
%         disp(' has already been renamed!')
    end
end


%% Copy to temp folder (due to access issues)
% Rename (update naming of ICA_clean file with sub-IDs)
destinationDir = '/data/MEG/Research/ECP/Semantic_Decision/temp_MEGnet'; % Directory with write access
newnames1 = cell(size(datafile_fif)); % Preallocate for efficiency
newnames2 = cell(size(datafile_fif)); % Preallocate for efficiency

for i = 1:length(datafile_fif)
    tkz = tokenize(datafile_fif{i}, '/');
    [pathstr, name, ext] = fileparts(datafile_fif{i});
    
    % Construct the new file name with destination directory
    newnames1{i} = fullfile(pathstr, [tkz{end-1}, '_', name, ext]);
    newnames2{i} = fullfile(destinationDir, [tkz{end-1}, '_', name, ext]);
        newnames{i} = fullfile(pathstr, [tkz{end-1}(1:end), '_', tkz{end}]);


    % Check if the file has already been copied to the new location
    if ~exist(newnames1{i}, 'file') && ~exist(newnames2{i}, 'file')
        disp(['Copying and renaming: ', datafile_fif{i}])
        commandStr = ['cp ', datafile_fif{i}, ' ', newnames2{i}];
%         pause,
        [status, cmdout] = system(commandStr);
        if status == 0
            disp('File copied and renamed.');
        else
            disp(['Error in copying file: ', cmdout]);
        end
    else
        disp([newnames2{i}, ' has already been copied and renamed!'])
    end
end


%%
% clc
% %- listing raw data
% tag = 'SD'; %- Semantic decision
% d = rdir([datadir,['/*/ICA/*/ica_clean.fif']]);
% 
% clear datafolder datafile sub subj_all
% datafile_fif = [];
% for i=1:length(d)
%     [pathstr, name] = fileparts(d(i).name);
%     datafolder{i} = pathstr;
%     datafile{i} = d(i).name;
%     Index = strfind(datafile{i}, '/');
%     sub{i} = datafile{i}(Index(7)+1:Index(8)-1);
%     subj_all{i} = [num2str(i), ': ', sub{i}];
% end
% datafile_fif = vertcat(datafile_fif,datafile);
% disp(datafile_fif')
% 
% %- Rename (update naeming of ICA_clean file with sub-IDs)
% newnames = [];
% for i=1:length(datafile_fif)
%     tkz = tokenize(datafile_fif{i},'/');
%     [pathstr, name] = fileparts(datafile_fif{i});
%     newnames{i} = fullfile(pathstr, [tkz{end-1}(1:end), '_', tkz{end}]);
%     cd(pathstr)
%     if isempty(rdir(newnames{i}))
%         copyfile(datafile_fif{i},newnames{i},'f');
%         disp(datafile_fif{i})
%         disp('copied.')
%     else
%         disp(datafile_fif{i})
%         disp(' has already been renamed!')
%     end
% end
% 
% %% Copy and rename the ICA_clean data
% clc
% %- listing raw data
% tag = 'SD'; %- Semantic decision
% d = rdir([datadir,['/*/ICA/*/ica_clean.fif']]);
% 
% clear datafolder datafile sub subj_all
% datafile_fif = [];
% for i=1:length(d)
%     [pathstr, name] = fileparts(d(i).name);
%     datafolder{i} = pathstr;
%     datafile{i} = d(i).name;
%     Index = strfind(datafile{i}, '/');
%     sub{i} = datafile{i}(Index(7)+1:Index(8)-1);
%     subj_all{i} = [num2str(i), ': ', sub{i}];
% end
% datafile_fif = vertcat(datafile_fif,datafile);
% disp(datafile_fif')
% 
% %- Rename (update naeming of ICA_clean file with sub-IDs)
% newnames = [];
% for i=1:length(datafile_fif)
%     tkz = tokenize(datafile_fif{i},'/');
%     [pathstr, name] = fileparts(datafile_fif{i});
%     newnames{i} = fullfile(pathstr, [tkz{end-1}(1:end), '_', tkz{end}]);
%     cd(pathstr)
%     if isempty(rdir(newnames{i}))
%         disp(datafile_fif{i})
%         pause,
%         copyfile(datafile_fif{i},newnames{i},'f');
%         disp('copied.')
%     else
%         disp(datafile_fif{i})
%         disp(' has already been renamed!')
%     end
% end




