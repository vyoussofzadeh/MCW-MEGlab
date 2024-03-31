
%% Test run for MEGnet

%% Initialize and set data directory
datadir = '/group/jbinder/ECP/MEG/MEG_Work';

% Search for raw .fif files in the specified directory structure
rawfiflist = rdir(fullfile(datadir,'/*/tSSS/*SD*raw.fif'));

% Clear Command Window for a clean output
clc

% Loop through each .fif file found
for i = 1:length(rawfiflist)
    % Extract the full path of the current .fif file
    rawfiflist_name = rawfiflist(i).name;
    % Extract the directory path and file name (without extension)
    [pathstr, name, ~] = fileparts(rawfiflist_name);
    
    % Determine the parent directory for the current file
    parentDir = fileparts(pathstr);
    
    % Find the occurrence of '_run' or '_Run' to identify run number
    runIdentifier = regexp(name, '_run|_Run', 'ignorecase', 'once');
    if isempty(runIdentifier)
        continue; % Skip if no run identifier is found
    end
    run_num = name(runIdentifier+1:runIdentifier+4);
    
    % Create 'ICA' directory if it does not exist
    icaDir = fullfile(parentDir, 'ICA');
    if ~exist(icaDir, 'dir'), mkdir(icaDir), end
    
    % Check if ICA cleaned data already exists to avoid re-processing
    if isempty(rdir(fullfile(icaDir, ['*', run_num, '*/ica_clean.fif'])))
        % Define the command for processing with MEGnet
        command = sprintf('ICA.py -filename %s/tSSS/%s -results_dir %s -line_freq 60', parentDir, [name, '.fif'], icaDir);
        disp(['MEGnet Processing of ', name])
        % Execute the command
        system(command);
    else
        disp(['ICA data for ', name, ' exists, skipped!']);
    end
end