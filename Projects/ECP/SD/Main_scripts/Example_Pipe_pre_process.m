%% The CID MEG

% MEG (pre)-processing pipeline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 08/30/2024

clear; clc, close('all'); warning off,

addpath('/group/bgross/work/CIDMEG/Analysis/CIDPipelines/functions/External')
addpath('/opt/mne_matlab/matlab')
datadir = '/group/bgross/work/CIDMEG/ECOG_MEG_data/MEG';
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/BS_additions/BS_event_read')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new/')
addpath('/group/bgross/work/CIDMEG/Analysis/CIDPipelines/functions')

ft_path = '/opt/matlab_toolboxes/ft_packages/fieldtrip_latest';
addpath(ft_path);
ft_defaults

%% MEGnet-ICA analysis - Part 1
cd(datadir)
rawfiflist = rdir(fullfile(datadir,'/mcw*/*/tsss/Block*raw.fif'));

clc
for i=1:length(rawfiflist)
    rawfiflist_name = rawfiflist(i).name;
    idx = strfind(rawfiflist_name,'/');
    [pathstr, name, ext] = fileparts(rawfiflist_name);
    data_filename = rawfiflist_name(idx(end)+1:end);
    cd(pathstr),
    
    idx1 = strfind(data_filename,'Block');
    idx2 = strfind(data_filename,'_');
    run_num = data_filename(idx1+5:idx2-1);
    if exist('ICA', 'file') == 0, mkdir('ICA'), end
    org_dir_name = ['./ICA/',data_filename(1:end-4)];
    new_dir_name = fullfile(pathstr, 'ICA', sprintf('Run_%s', run_num));
    
    if ~exist(new_dir_name, 'file')
        command = ['/opt/miniconda3/envs/megnet/bin/python /opt/miniconda3/bin/ICA.py -filename $(pwd)/', data_filename, ' -results_dir $(pwd)/ICA -line_freq  60'];
        disp(['MEGnet Processing of ', data_filename])
        system(command)
    else
        disp(['ICA data for ', name,' exist, skipped!'])
    end
    
    % Directory renaming logic
    if exist(org_dir_name, 'dir') && exist(fullfile(org_dir_name, 'ica_clean.fif'), 'file') && ~exist(new_dir_name, 'dir')
        movefile(org_dir_name, new_dir_name);
        disp(['ICA directory ', data_filename(1:end-4), ' renamed to Run_', run_num]);
    end
end

%% MEGnet-ICA analysis - Part 2 data
cd(datadir)
rawfiflist = rdir(fullfile(datadir,'/mcw*/*/tsss/Run*cidMEG*raw.fif'));

clc
for i=1:length(rawfiflist)
    
    rawfiflist_name = rawfiflist(i).name;
    idx = strfind(rawfiflist_name,'/');
    [pathstr, name, ext] = fileparts(rawfiflist_name);
    data_filename = rawfiflist_name(idx(end)+1:end);
    cd(pathstr),
    %     cd ..
    idx1 = strfind(data_filename,'Run');
    idx2 = strfind(data_filename,'_');
    run_num = data_filename(idx1+3:idx2-1);
    if exist('ICA', 'file') == 0, mkdir('ICA'), end
    org_dir_name = ['./ICA/',data_filename(1:end-4)];
    new_dir_name = fullfile(pathstr, 'ICA', sprintf('Run_%s', run_num));
    if ~exist(new_dir_name, 'file')
        command = ['/opt/miniconda3/envs/megnet/bin/python /opt/miniconda3/bin/ICA.py -filename $(pwd)/', data_filename, ' -results_dir $(pwd)/ICA -line_freq  60'];
        disp(['MEGnet Processing of ', data_filename])
        system(command)
    else
        disp(['ICA data for ', name,' exist, skipped!'])
    end
    
    % Directory renaming logic
    if exist(org_dir_name, 'dir') && exist(fullfile(org_dir_name, 'ica_clean.fif'), 'file') && ~exist(new_dir_name, 'dir')
        movefile(org_dir_name, new_dir_name);
        disp(['ICA directory ', data_filename(1:end-4), ' renamed to Run_', run_num]);
    end
end

%% Read raw data
% clc
cd(datadir)
tag = 'CID';

run = [];
savefile = ['dataprocesslog_',tag,'.mat'];

clc
clear datafolder datafile subj_all sub
datafile_fif = [];
d = rdir([datadir,['/**/ICA/*','*/ica_clean.fif']]);
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafile{i} = d(i).name;
    
    parts = strsplit(pathstr, '/');
    subjectID = parts{8};  % Adjust the index based on your path structure
    
    Datarun = parts{end};
    subj_all{i} = [num2str(i), ': ', subjectID,'_', Datarun];
end
datafile_fif = vertcat(datafile_fif,datafile);
datafile_fif = datafile_fif';

disp(datafile_fif)
disp(subj_all')

%% BS
% bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3';
bs_path = '/data/MEG/Vahab/Github/brainstorm3';

addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

BS_dir = '/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/Brainstorm_db/CID';
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

%% Import raw data into BS
clc
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/BS_additions/BS_event_read/cid projects')

L_data =length(datafile_fif);
not_imported = [];
for ii = 1:L_data
    
    [pathname, name] = fileparts(datafile_fif{ii});
    [pathstr2, name2] = fileparts(pathstr);
    
    pattern = '\w+_v\d+';
    extractedStr = regexp(subj_all{ii}, pattern, 'match');
    extractedStr = extractedStr{1};
    
    datadir_sub = fullfile(BS_dir,'data/',extractedStr);
    datadir_sub = strrep(datadir_sub, ' ', '');
    
    tkz = tokenize(pathname,'/');
    sub_sel = tkz{8};
    cd(datadir_sub);
    
    idx = strfind(datadir_sub,'/');
    
    % Extracting the Run Number
    runNumberPattern = 'Run_(\d+)';
    [~, ~, ~, match] = regexp(pathname, runNumberPattern);
    runNumber = match{1}; % Convert the string to a numeric value
    
    if  isempty(dir(fullfile(['./@rawica_',runNumber, '_clean'],'/channel_*.mat'))) && ...
            isempty(dir(fullfile(['./@raw',runNumber],'/channel_*.mat')))
        
        iSubject = find(contains(unq_bs_subj, sub_sel)==1);
        RawFiles = datafile_fif{ii};
        disp(RawFiles)
        %         pause,
        
        [fPath, ~] = bst_fileparts(RawFiles);[fPath, fBase] = bst_fileparts(fPath);
        cid_import_raw(RawFiles, iSubject, fBase); % no GUI
    end
end

%% Preprocess raw data
d = rdir(fullfile(BS_data_dir,'/mcw*/@raw*/*raw_ica_clean.mat'));

for i=1:length(d)
    
    [pathstr, name] = fileparts(d(i).name);
    [pathstr2, name2] = fileparts(pathstr);
    [pathstr3, sub_sel] = fileparts(pathstr2);
    sFiles = {fullfile(sub_sel, name2, [name, '.mat'])};
    
    pathStr = d(i).name;
    
    % Extracting the Subject ID
    subjectIDPattern = 'mcwa\d+_v\d+';
    [~, ~, ~, match] = regexp(pathStr, subjectIDPattern);
    subjectID = match{1};
    
    % Extracting the Run Number
    runNumberPattern = 'Run_(\d+)';
    [~, ~, ~, match] = regexp(pathStr, runNumberPattern);
    runNumber = match{1}; % Convert the string to a numeric value
    
    % Display the extracted information
    fprintf('Subject ID: %s\n', subjectID);
    fprintf('Run Number: %s\n', runNumber);
    
    newFolderName = fullfile(pathstr2, ['@rawica_', runNumber, '_clean']);
    
    if ~exist(newFolderName, 'dir')
        
        db_reload_database('current',1);
        pause(3);
        cd(pathstr2)
        
        disp('preprocessing,')
        disp(sFiles)
        sFiles_out = cid_bs_preprocess(sFiles);
        pause(4)
        
        redundantfolder = fullfile(pathstr2, [sFiles_out.Condition]);
        if exist(redundantfolder, 'dir')
            rmdir(redundantfolder,'s')
        end
        
        oldFolderName = fullfile(pathstr2, [sFiles_out.Condition, '_clean']);
        if exist(oldFolderName, 'dir')
            movefile(oldFolderName, newFolderName);
        end
    end
end

disp('Preprocessing was completed.')

%% Import epoched data
db_reload_database('current',1)
cd(BS_data_dir)
d = rdir(fullfile(BS_data_dir,'/mcw*/@raw*/*band_clean.mat'));

clc
for i=1:length(d)
    
    [pathstr, name] = fileparts(d(i).name);
    [pathstr2, name2] = fileparts(pathstr);
    [pathstr3, sub_sel] = fileparts(pathstr2);
    sFiles = {fullfile(sub_sel, name2, [name, '.mat'])};
    
    pathStr = d(i).name;
    
    % Extracting the Subject ID
    subjectIDPattern = 'mcwa\d+_v\d+';
    [~, ~, ~, match] = regexp(pathStr, subjectIDPattern);
    subjectID = match{1};
    
    % Extracting the Run Number
    runNumberPattern = 'Run_(\d+)';
    [~, ~, ~, match] = regexp(pathStr, runNumberPattern);
    runNumber = match{1}; % Convert the string to a numeric value
    
    % Display the extracted information
    fprintf('Subject ID: %s\n', subjectID);
    fprintf('Run Number: %s\n', runNumber);
    
    cd(pathstr2)
    newFolderName = fullfile(pathstr2, [subjectID, '_', runNumber]);
    
    if ~exist(newFolderName, 'dir')
        
        disp(sFiles)
        % Define epoching parameters
        %         EpochTime = [-0.3, 0.6];  % Epoch from -100 ms to +300 ms around the event
        EpochTime = [0, 1];  % Epoch from -100 ms to +300 ms around the event
        
        subjectName = runNumber;
        
        db_reload_database('current',1);
        
        switch subjectID
            case 'mcwa038_v1'
                eventname = 'STI102';
            case 'mcwa065_v1'
                eventname = 'STI101';
        end
        
        % Create epochs
        sFiles_out = bst_process('CallProcess', 'process_import_data_event', sFiles{1}, [], ...
            'subjectname', subjectID, ...
            'condition', '', ...
            'eventname', eventname, ...
            'timewindow', [], ...
            'epochtime', EpochTime, ...
            'createcond', 0, ...
            'ignoreshort', 1, ...
            'usectfcomp', 1, ...
            'usessp', 1, ...
            'freq', 1000, ...
            'baseline',    [-0.3, -0.001]);
        
        pause(3);
        cd(pathstr2)
        
        oldFolderName = fullfile(pathstr2, [sFiles_out(1).Condition]);
        if exist(oldFolderName, 'dir')
            movefile(oldFolderName, newFolderName);
        end
    end
end

%% Est. head model
d = rdir(fullfile(BS_data_dir,'/mcw*/mcwa*/channel_vectorview306_acc1.mat'));

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

%% Source Analysis - LCMV
d = rdir(fullfile(BS_data_dir,'/mcw*/mcwa*/channel_vectorview306_acc1.mat'));

for i=1:length(d)
    
    [pathstr, name] = fileparts(d(i).name);
    dd_results = rdir(fullfile(pathstr,'results*.mat'));
    results = [];
    for j=1:length(dd_results)
        cmt = load(dd_results(j).name);
        results{j} = cmt.Comment;
    end
    
    cd(pathstr)
    
    if contains(pathstr, 'mcwa038_v1')
        dtag = 'STI102';
    elseif contains(pathstr, 'mcwa065_v1')
        dtag = 'STI101';
    end
    
    dd = rdir(fullfile(pathstr,['data_', dtag, '*_trial*.mat']));
    
    if isempty(results)
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
        cfg = []; cfg.comment = 'source'; bs_lcmv(cfg, sFiles)
    end
end

%% Project to Default Template
destSurfFile = '@default_subject/tess_cortex_pial_low.mat';
% Processing for projecting to default template
load(protocol);
Subj = ProtocolSubjects.Subject;
ProtocolSubjects = []; k=1;
for i=1:length(Subj)
    if contains(Subj(i).Name, 'mcw')
        ProtocolSubjects{k} =  Subj(i).Name; k=1+k;
    end
end

subj = ProtocolSubjects;

etag = 'mcwa';

for ii = 1:length(subj)
    
    switch subj{ii}
        case 'mcwa038_v1'
            eventname = 'STI102';
        case 'mcwa065_v1'
            eventname = 'STI101';
    end
    cd(fullfile(BS_data_dir,subj{ii}))
    dd = rdir(['./',etag,'*/results_PNAI_*.mat']);
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
        d = rdir(['./',fileparts(sFiles1{j}), '/data_', eventname, '*_average*.mat']);
        disp(sFiles1{j})
        if length(d) > 1
            for jj=1:length(d)
                tmp2 = load(d(jj).name);
                disp([num2str(jj), ': ', tmp2.Comment])
            end
            del_sel = input('sel_del:');
            delete(d(del_sel).name)
            d = rdir(['./',fileparts(sFiles1{j}), '/data_', eventname, '*_average*.mat']);
        end
        sFiles2 = ['link|',sFiles1{j}, '|',d.name(3:end)];
        
        tkz = tokenize(sFiles1{j},'/');
        d_name = [];
        dd = rdir(fullfile(BS_data_dir, '/Group_analysis', tkz{2}, 'results_*.mat'));
        if isempty(dd)
            sFiles = bst_project_sources({sFiles2}, destSurfFile, 0, 1);
        end
    end
end

%% Extract ROI sources, Brainnetome
disp('Processing source files from Group_analysis to derive time series...');

% Define the Group_analysis directory
group_analysis_dir = fullfile(BS_data_dir, 'Group_analysis');
results_files = rdir(fullfile(group_analysis_dir, '**/results_PNAI_*.mat'));

if isempty(results_files)
    error('No results_PNAI_*.mat files found in Group_analysis directory.');
end

% Define the save directory and ensure it exists
savedir = '/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/ROI_timeseries';
if ~exist(savedir, 'dir')
    mkdir(savedir);
    disp(['Created directory: ', savedir]);
end

% Loop through each results file
for iFile = 1:length(results_files)
    results_file = results_files(iFile).name;
    
    % Process: Scout time series: [72 scouts]
    sFiles = bst_process('CallProcess', 'process_extract_scout', {results_file}, [], ...
        'timewindow',     [0, 1], ...
        'scouts',         {'Brainnetome_association', {'A20cl_L', 'A20cl_R', 'A20cv_L', 'A20cv_R', 'A20il_L', 'A20il_R', 'A20iv_L', 'A20iv_R', 'A20r_L', 'A20r_R', 'A20rv_L', 'A20rv_R', 'A21c_L', 'A21c_R', 'A21r_L', 'A21r_R', 'A22c_L', 'A22c_R', 'A22r_L', 'A22r_R', 'A37dl_L', 'A37dl_R', 'A37elv_L', 'A37elv_R', 'A37lv_L', 'A37lv_R', 'A37mv_L', 'A37mv_R', 'A37vl_L', 'A37vl_R', 'A38l_L', 'A38l_R', 'A38m_L', 'A38m_R', 'A41/42_L', 'A41/42_R', 'TE1.0/TE1.2_L', 'TE1.0/TE1.2_R', 'TI_L', 'TI_R', 'aSTS_L', 'aSTS_R', 'cpSTS_L', 'cpSTS_R', 'rpSTS_L', 'rpSTS_R', 'A7c_L', 'A7c_R', 'A7ip_L', 'A7ip_R', 'A7m_L', 'A7m_R', 'A7pc_L', 'A7pc_R', 'A7r_L', 'A7r_R', 'A39c_L', 'A39c_R', 'A39rd_L', 'A39rd_R', 'A39rv_L', 'A39rv_R', 'A40c_L', 'A40c_R', 'A40rd_L', 'A40rd_R', 'A40rv_L', 'A40rv_R', 'A5l_L', 'A5l_R', 'A5m_L', 'A5m_R'}}, ...
        'flatten',        1, ...
        'scoutfunc',      'mean', ...  % Mean
        'pcaedit',        struct(...
        'Method',         'pca', ...
        'Baseline',       [-0.3, 2], ...
        'DataTimeWindow', [-0.3, 2], ...
        'RemoveDcOffset', 'file'), ...
        'isflip',         1, ...
        'isnorm',         0, ...
        'concatenate',    0, ...
        'save',           1, ...
        'addrowcomment',  1, ...
        'addfilecomment', []);
    
    
    disp(['Processing file: ', results_file]);
    
    % Load the results file
    sResults = in_bst_results(results_file, 1);
    
    % Check if source data exists in the results file
    if ~isfield(sResults, 'ImageGridAmp') || isempty(sResults.ImageGridAmp)
        disp(['No source data in file: ', results_file]);
        continue;
    end
    
    % Extract identifiers from HeadModelFile
    headModelPath = sResults.HeadModelFile;  % e.g., 'mcwa065_v1/mcwa065_v1_Run_9/headmodel_surf_os_meg.mat'
    tokens = regexp(headModelPath, '(?<subjID>[^/]+)/.*(?<runID>Run_\d+)', 'names');
    
    % Use tokens for naming
    if ~isempty(tokens)
        subjID = tokens.subjID;   % e.g., 'mcwa065_v1'
        runID = tokens.runID;     % e.g., 'Run_9'
    else
        warning(['Unable to parse HeadModelFile: ', headModelPath]);
        subjID = 'UnknownSubject';
        runID = 'UnknownRun';
    end
    
    % Load the Brainnetome atlas for scout mapping
    [sScout, ProtocolInfo] = do_convertBrainnetomeScout(sResults);
    if isempty(sScout)
        disp(['No valid Brainnetome scouts found for file: ', results_file]);
        continue;
    end
    
    % Derive time series for all scouts
    %     scout_time_series = cid_process_extract_scout_time_series(sResults, sScout, 'mean');
    vs = load(fullfile(BS_data_dir, sFiles.FileName));
    
    % Define filename for the subject-run
    output_filename = sprintf('%s_%s_ROI_timeseries.mat', subjID, runID);
    output_filepath = fullfile(savedir, output_filename);
    
    % Save all ROIs for the current subject-run into one file
    timepoints = sResults.Time;  % Add time vector
    save(output_filepath, '-struct','vs');
    disp(['Saved ROI time series for: ', subjID, ' ', runID, ' to ', output_filepath]);
end

disp('All ROI time series have been processed and saved.');

cd(savedir)

%% Process Source Files from Group_analysis - all verticies, optional
if extractallverticies == 1
    disp('Processing source files from Group_analysis to derive time series...');
    
    % Define the Group_analysis directory
    group_analysis_dir = fullfile(BS_data_dir, 'Group_analysis');
    results_files = rdir(fullfile(group_analysis_dir, '**/results_PNAI_*.mat'));
    
    if isempty(results_files)
        error('No results_PNAI_*.mat files found in Group_analysis directory.');
    end
    
    % Define the save directory and ensure it exists
    savedir = '/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/ROI_timeseries';
    if ~exist(savedir, 'dir')
        mkdir(savedir);
        disp(['Created directory: ', savedir]);
    end
    
    % Loop through each results file
    for iFile = 1:length(results_files)
        results_file = results_files(iFile).name;
        disp(['Processing file: ', results_file]);
        
        % Load the results file
        sResults = in_bst_results(results_file, 1);
        
        % Check if source data exists in the results file
        if ~isfield(sResults, 'ImageGridAmp') || isempty(sResults.ImageGridAmp)
            disp(['No source data in file: ', results_file]);
            continue;
        end
        
        % Extract identifiers from HeadModelFile
        headModelPath = sResults.HeadModelFile;  % e.g., 'mcwa065_v1/mcwa065_v1_Run_9/headmodel_surf_os_meg.mat'
        tokens = regexp(headModelPath, '(?<subjID>[^/]+)/.*(?<runID>Run_\d+)', 'names');
        
        % Use tokens for naming
        if ~isempty(tokens)
            subjID = tokens.subjID;   % e.g., 'mcwa065_v1'
            runID = tokens.runID;     % e.g., 'Run_9'
        else
            warning(['Unable to parse HeadModelFile: ', headModelPath]);
            subjID = 'UnknownSubject';
            runID = 'UnknownRun';
        end
        
        % Load the Brainnetome atlas for scout mapping
        [sScout, ProtocolInfo] = do_convertBrainnetomeScout(sResults);
        if isempty(sScout)
            disp(['No valid Brainnetome scouts found for file: ', results_file]);
            continue;
        end
        
        % Derive time series for all scouts
        scout_time_series = cid_process_extract_scout_time_series(sResults, sScout, 'all');
        
        % Define filename for the subject-run
        output_filename = sprintf('%s_%s_ROI_timeseries.mat', subjID, runID);
        output_filepath = fullfile(savedir, output_filename);
        
        % Save all ROIs for the current subject-run into one file
        timepoints = sResults.Time;  % Add time vector
        save(output_filepath, 'scout_time_series', 'timepoints');
        disp(['Saved ROI time series for: ', subjID, ' ', runID, ' to ', output_filepath]);
    end
    
    disp('All ROI time series have been processed and saved.');
    
end
