
function cleanup_MaxFilter_tSSS_forWrapper_ver2(subjIDpath, outDirpath)
% subjIDpath: folder that contains *raw.fif
% outDirpath: folder where outputs should be written (tSSS outputs)

% ---- normalize inputs (R2018a safe)
subjIDpath = char(subjIDpath);
outDirpath = char(outDirpath);

fprintf('Input file path: %s\n', subjIDpath);
fprintf('Output files will be saved in: %s\n', outDirpath);

if exist(subjIDpath,'dir') ~= 7
    error('Input folder not found: %s', subjIDpath);
end
if exist(outDirpath,'dir') ~= 7
    mkdir(outDirpath);
end

% --- Create meg_clinic XML templates (session.xml, raw.xml)
make_megclinic_xml_templates(subjIDpath, outDirpath);

do_analysis.ecg = 1;
do_analysis.eog = 0;
do_analysis.ongoing = 1;

% ---- Add paths
pathList = { ...
    '/MEG_data/megclinic', ...
    '/usr/local/MATLAB_Tools/brainstorm3', ...
    '/usr/local/MNE-2.7.0-3106-Linux-x86_64/share/matlab', ...
    '/usr/local/MATLAB_Tools/mne', ...
    '/neuro/bin' ...
};
for ii = 1:numel(pathList)
    if exist(pathList{ii}, 'dir') == 7
        addpath(genpath(pathList{ii}));
    else
        warning('Path not found: %s', pathList{ii});
    end
end

% ---- Work in input folder
cd(subjIDpath);
files = dir('*raw.fif');
if isempty(files)
    fprintf('No *raw.fif found in %s\n', subjIDpath);
    return;
end

for j = 1:numel(files)
    fileName = files(j).name;
    baseFileName = strrep(fileName, '.fif', '');

    % Naming structs:
    % IMPORTANT: baseDir here should be the OUTPUT destination for tSSS/cleanup products
    FileNames_tsss = generateFileNames_tsss_cleanup(outDirpath, fileName(1:end-8));
    FileNames_sss  = generateFileNames_sss_nocleanup(outDirpath, fileName(1:end-8));

    isEmpty = ~isempty(regexpi(baseFileName, 'emptyroom'));

    if isEmpty
        % emptyroom: do SSS only
        if ~exist(fullfile(FileNames_sss.filelocation, FileNames_sss.cleanFileName),'file')
            disp(['Empty-room detected. Running SSS for: ', baseFileName]);
            runMaxFilterSSS(subjIDpath, outDirpath, files(j));  % PASS outDirpath
            disp('Skipping SSP cleanup for empty-room.');
        else
            disp(['Skipping ', baseFileName, ' (SSS exists).']);
        end
        continue;
    end

    % non-emptyroom: do tSSS then cleanup
    if ~exist(fullfile(FileNames_tsss.filelocation, FileNames_tsss.cleanFileName),'file')
        disp(['Preprocessing: ', baseFileName]);

        if exist(FileNames_tsss.filename,'file') ~= 2
            runMaxFilterTSSS(subjIDpath, outDirpath, files(j)); % PASS outDirpath
        else
            fprintf('✓ %s — tSSS already exists, skipping MaxFilter.\n', baseFileName);
        end

        if do_analysis.ongoing, nogui_remove_ongoing_artifact(FileNames_tsss); end
        if do_analysis.ecg,     nogui_remove_ecg_artifact(FileNames_tsss); end
        if do_analysis.eog,     nogui_remove_eog_artifact(FileNames_tsss); end

        % apply projections to make clean file
        cleanFile = fullfile(FileNames_tsss.filelocation, FileNames_tsss.cleanFileName);
        if ~exist(cleanFile,'file')
            appliedProj = {};

            if do_analysis.ecg
                if exist(fullfile(FileNames_tsss.filelocation, FileNames_tsss.ecgProjFileName),'file')
                    appliedProj{end+1} = FileNames_tsss.ecgProjFileName;
                end
            end
            if do_analysis.eog
                if exist(fullfile(FileNames_tsss.filelocation, FileNames_tsss.eogProjFileName),'file')
                    appliedProj{end+1} = FileNames_tsss.eogProjFileName;
                end
            end

            if ~isempty(appliedProj)
                projString = '';
                for k = 1:numel(appliedProj)
                    projString = [projString ' --proj ' appliedProj{k}]; %#ok<AGROW>
                end

                cmd = sprintf('mne_process_raw --cd %s --raw %s %s --projon --save %s --filteroff', ...
                    FileNames_tsss.filelocation, FileNames_tsss.filename, projString, FileNames_tsss.cleanFileName);

                [status, cmdout] = unix(cmd);
                if status
                    error('Error processing MNE command:\n%s', cmdout);
                end
            end
            disp('Completed.');
        end
    else
        disp(['Skipping ', baseFileName, ' (already cleaned).']);
    end
end
end

function runMaxFilterTSSS(inDir, outDir, fileStruct)
fileName = fileStruct.name;
inputFile = fullfile(inDir, fileName);
subDirName = fileName(1:end-8);

outputSubDir = fullfile(outDir, 'sss', subDirName);
if exist(outputSubDir,'dir') ~= 7
    mkdir(outputSubDir);
end
outputFile = fullfile(outputSubDir, [subDirName, '_raw_t_sss.fif']);

cmd = sprintf(['/neuro/bin/util/maxfilter -gui -f %s -o %s ' ...
    '-ctc /neuro/databases/ctc/ct_sparse.fif ' ...
    '-cal /neuro/databases/sss/sss_cal.dat ' ...
    '-autobad off -st 10 -corr 0.9 -force'], inputFile, outputFile);

[status, cmdout] = system(cmd);
if status ~= 0
    fprintf('Error running MaxFilter for %s:\n%s\n', fileName, cmdout);
else
    fprintf('Successfully processed %s\n', fileName);
end
end

function runMaxFilterSSS(inDir, outDir, fileStruct)
fileName = fileStruct.name;
inputFile = fullfile(inDir, fileName);
subDirName = fileName(1:end-8);

outputSubDir = fullfile(outDir, 'sss', subDirName);
if exist(outputSubDir,'dir') ~= 7
    mkdir(outputSubDir);
end
outputFile = fullfile(outputSubDir, [subDirName, '_raw_sss.fif']);

cmd = sprintf(['/neuro/bin/util/maxfilter -gui -f %s -o %s ' ...
    '-ctc /neuro/databases/ctc/ct_sparse.fif ' ...
    '-cal /neuro/databases/sss/sss_cal.dat ' ...
    '-autobad off -force'], inputFile, outputFile);

[status, cmdout] = system(cmd);
if status ~= 0
    fprintf('Error running MaxFilter for %s:\n%s\n', fileName, cmdout);
else
    fprintf('Successfully processed %s\n', fileName);
end
end

function FileNames = generateFileNames_tsss_cleanup(baseDir, baseFileName)

filelocation = fullfile(baseDir, 'sss', baseFileName);

% Initialize the FileNames structure with paths
FileNames = struct();
FileNames.filelocation = filelocation;
FileNames.filename = fullfile(filelocation, [baseFileName, '_raw_t_sss.fif']);
FileNames.ogProjFileName = fullfile(filelocation, [baseFileName, '_raw_t_sss_ongoing-proj.fif']);
FileNames.ogCleanFileName = [baseFileName, '_raw_t_sss_ongoingClean_raw.fif'];

FileNames.eventFileName = [baseFileName, '_raw_t_sss_ecg-eve.fif'];
FileNames.ecgProjFileName =  [baseFileName, '_raw_t_sss_ecg-proj.fif'];

FileNames.eogEventFileName = [baseFileName, '_raw_t_sss_eog-eve.fif'];
FileNames.eogProjFileName =  [baseFileName, '_raw_t_sss_eog-proj.fif'];

FileNames.cleanFileName = [baseFileName, '_raw_t_sss_ecgClean_raw.fif'];

FileNames.ecgCleanFileName = [baseFileName, '_raw_t_sss_ecgClean_raw.fif'];
FileNames.eogCleanFileName =  [baseFileName, '_raw_t_sss_eogClean_raw.fif'];

end

function FileNames = generateFileNames_sss_nocleanup(baseDir, baseFileName)

filelocation = fullfile(baseDir, 'sss', baseFileName);

% Initialize the FileNames structure with paths
FileNames = struct();
FileNames.filelocation = filelocation;
FileNames.filename = fullfile(filelocation, [baseFileName, '_raw_sss.fif']);
FileNames.cleanFileName = [baseFileName, '_raw_sss.fif'];

end

function nogui_make_proj_operator(directory, raw, eve, type)
% make_proj_operator: create ssp vectors and save to a file
%
% USAGE:    make_proj_operator(directory, raw, eve, type)   for making projections from an event file
%           make_proj_operator(directory, raw, [], type)    for making projections over entire file
%
% INPUT:    directory = location of the files
%           raw = name of the file to use for creating projections
%           eve = name of the event file
%           type = [ECG, ONGOING] type of projections to create
%
%
% Author: Elizabeth Bock, 2009
% --------------------------- Script History ------------------------------
% EB 27-AUG-2009  Creation
% -------------------------------------------------------------------------
% logFile = GUI.MCLogFile;
% config = ArtifactClean.CleanConfig;

switch (type)
    case 'ECG'
        if strfind(raw, '_raw_sss')
            projTag = strcat('_raw_sss',char('_ecg-proj'));
        else
            projTag = char('_ecg-proj');
        end
        start =  '-0.08'; %char(config.ECG_TMIN);
        stop = '0.08'; % char(config.ECG_TMAX);
        nMag = '1'; %char(config.ECG_NMAG);
        nGrad = '1'; %char(config.ECG_NGRAD);
        filterOn = true(1); %config.INCLUDE_ECG_FILTERING;
        
        if filterOn
            hpf = '10'; %char(config.ECG_HPFILTER);
            lpf = '40'; %char(config.ECG_LPFILTER);
            command = ['mne_process_raw --cd ' directory ' --raw ' raw ' --events ' eve ' --makeproj --projtmin ' start ' --projtmax ' stop ' --saveprojtag ' projTag ' --projnmag ' nMag ' --projngrad ' nGrad ' --projevent 999 --highpass ' hpf ' --lowpass ' lpf ' --digtrigmask 0'];
            %             logFile.write(['command: ' command]);
            [status,w] = unix(command);
            %             logFile.write(w);
        else
            command = ['mne_process_raw --cd ' directory ' --raw ' raw ' --events ' eve ' --makeproj --projtmin ' start ' --projtmax ' stop ' --saveprojtag ' projTag ' --projnmag ' nMag ' --projngrad ' nGrad ' --projevent 999 --filteroff --digtrigmask 0'];
            %             logFile.write(['command: ' command]);
            [status,w] = unix(command);
            %             logFile.write(w);
        end
        
        %  For names containing DefaultHead_sss
        if strfind(raw, '_defaultHead_sss.fif')
            projFileName = strrep(raw, '_defaultHead_sss.fif', '_defaultHead_raw_sss_ecg-proj.fif');
            if exist(projFileName, 'file')
                newProjFileName = strrep(projFileName, '_defaultHead_raw_sss_ecg-proj.fif', '_defaultHead_sss_ecg-proj.fif');
                movefile(projFileName, newProjFileName);
            end
        end
        
    case 'ONGOING'
        projTag = strcat('_raw_sss_ongoing-proj');
        nMag = '3';%char(config.ONGOING_NMAG);
        nGrad = '3';%char(config.ONGOING_NGRAD);
        filterOn = true(1); %config.INCLUDE_OG_FILTERING;
        
        if filterOn
            hpf = '1.5'; % char(config.OG_HPFILTER);
            lpf = '5'; %char(config.OG_LPFILTER);
            command = ['mne_process_raw --cd ' directory ' --raw ' raw ' --makeproj --saveprojtag ' projTag ' --projnmag ' nMag ' --projngrad ' nGrad ' --highpass ' hpf ' --lowpass ' lpf ' --digtrigmask 0 --projgradrej -1 --projmagrej -1'];
            [status,w] = unix(command);
        else
            command = ['mne_process_raw --cd ' directory ' --raw ' raw ' --makeproj --saveprojtag ' projTag ' --projnmag ' nMag ' --projngrad ' nGrad ' --filteroff --digtrigmask 0 -projgradrej -1 -projmagrej -1'];
            [status,w] = unix(command);
        end
        
    case 'EOG'
        if strfind(raw, '_raw_sss')
            projTag = strcat('_raw_sss',char('_eog-proj'));
        else
            projTag = char('_eog-proj');
        end
        start =  '-0.2'; %char(config.EOG_TMIN);
        stop = '0.2'; %char(config.EOG_TMAX);
        nMag = '1'; %char(config.EOG_NMAG);
        nGrad = '1'; %char(config.EOG_NGRAD);
        filterOn = true(1); %config.INCLUDE_EOG_FILTERING;
        
        events = mne_read_events(fullfile(directory,eve));
        eventNo = mode(double(events(:,3)));
        
        if filterOn
            hpf = '1.5'; % char(config.EOG_HPFILTER);
            lpf = '15'; % char(config.EOG_LPFILTER);
            command = ['mne_process_raw --cd ' directory ' --raw ' raw ' --events ' eve ' --makeproj --projtmin ' start ' --projtmax ' stop ' --saveprojtag ' projTag ' --projnmag ' nMag ' --projngrad ' nGrad ' --projevent ' num2str(eventNo) ' --highpass ' hpf ' --lowpass ' lpf ' --digtrigmask 0'];
            [status,w] = unix(command);
        else
            command = ['mne_process_raw --cd ' directory ' --raw ' raw ' --events ' eve ' --makeproj --projtmin ' start ' --projtmax ' stop ' --saveprojtag ' projTag ' --projnmag ' nMag ' --projngrad ' nGrad ' --projevent ' num2str(eventNo) ' --filteroff --digtrigmask 0'];
            [status,w] = unix(command);
        end
        
        %  For names containing DefaultHead_sss
        if strfind(raw, '_defaultHead_sss.fif')
            projFileName = strrep(raw, '_defaultHead_sss.fif', '_defaultHead_raw_sss_ecg-proj.fif');
            if exist(projFileName, 'file')
                newProjFileName = strrep(projFileName, '_defaultHead_raw_sss_ecg-proj.fif', '_defaultHead_sss_ecg-proj.fif');
                movefile(projFileName, newProjFileName);
            end
        end
end
end

function nogui_remove_ongoing_artifact(FileNames)
% remove_ongoing_artifact: create an SSP projection for ongoing artifacts
%
% USAGE:    remove_ongoing_artifact(mc, FileNames)
%
% INPUT:    mc = megclinic instance
%           FileNames = Filename structure
%
%
% Author: Elizabeth Bock, 2009
% --------------------------- Script History ------------------------------
% EB 27-AUG-2009  Creation
% -------------------------------------------------------------------------

% -------------- SSP ------------------------------------------------------
% check to see if proj file already exists
checkfile = fullfile(FileNames.filelocation, FileNames.ogProjFileName);
projExists = exist(checkfile,'file');
if (projExists == 0)
    % create new ssp operators
    nogui_make_proj_operator(FileNames.filelocation, FileNames.filename, [], 'ONGOING');
    
    % find all _ongoing-proj.fif files
    projFile = dir(fullfile(FileNames.filelocation, '*_ongoing-proj.fif'));
    % find the most recent
    [date,index] = max([projFile.datenum]);
    projFile = fullfile(FileNames.filelocation,projFile(index).name);
    
    % Rename this file if different than the convention
    %     if ~strcmp(projFile, fullfile(FileNames.filelocation, FileNames.ogProjFileName))
    %         movefile(projFile, fullfile(FileNames.filelocation, FileNames.ogProjFileName));
    %     end
else
    disp('ongoing SSP complete')
end
end

function nogui_remove_eog_artifact(FileNames)
% remove_eog_artifact: Get eye blinks and create SSP projection
%
% USAGE:    remove_eog_artifact(mc, FileNames)
%
% INPUT:    mc = megclinic instance
%           FileNames = FileNames structure
%
%
% Author: Elizabeth Bock, 2010
% --------------------------- Script History ------------------------------
% EB 18-NOV-2010  Creation
% -------------------------------------------------------------------------

%% Get eye blink events
eogEveFile = fullfile(FileNames.filelocation, FileNames.eogEventFileName);
fileExists = exist(eogEveFile,'file');
if fileExists == 0
    fiffsetup = fiff_setup_read_raw(FileNames.filename);
    % Get eog
    eog = get_eog(fiffsetup);
    Options.sampRate = double(fiffsetup.info.sfreq);
    Options.firstSamp = double(fiffsetup.first_samp);
    % Find events
    Options.percentThresh = 80; % Threshold
    Options.ampMin = 85e-6;     % minimum amplitude
    Options.maxCross = 3;       % max threshold crossings
    Options.minBeats = 5;      % min blinks for ssp
    Options.corrVal = 0.5;      % correlation cutoff for sorting
    eog_events = blinkDetect(eog, Options);
    
    if ~isempty(eog_events)
        mneEvent(1,1) = Options.firstSamp;
        mneEvent(1,2) = Options.firstSamp/Options.sampRate;
        mneEvent(1,3:4) = 0;
        
        saveEvents = [mneEvent;eog_events];
        eventList(:,1) = saveEvents(:,1);
        eventList(:,2) = saveEvents(:,3);
        eventList(:,3) = saveEvents(:,4);
        
        % Write to .fif file
        mne_write_events(eogEveFile,eventList);
    end
end

% Calculate SSP
%check for existing projection and event files
projExists = exist(fullfile(FileNames.filelocation, FileNames.eogProjFileName),'file');
eveExists = exist(fullfile(FileNames.filelocation, FileNames.eogEventFileName), 'file');
if ~projExists
    if eveExists
        % create new ssp operators
        %         mc.setMessage(GUI.Config.M_MAKE_SSP);
        make_proj_operator(FileNames.filelocation, FileNames.filename, FileNames.eogEventFileName, 'EOG');
        
        projExists = exist(fullfile(FileNames.filelocation, FileNames.eogProjFileName),'file');
        % If the proj file does not exist, then check to be sure MNE named
        % the file correctly.  If not, change the file name.
        if ~projExists
            % find all -proj files for eog
            projFile = dir(fullfile(FileNames.filelocation, ['*' char(ArtifactClean.CleanConfig.EOG_PROJ) '.fif']));
            if isempty(projFile)
                % No proj files for eog exist
                %                 GUI.ErrorMessage(GUI.ErrorMessage.GENERIC_WARNING, sprintf('%s\n%s','Eye blink projections not created.',' See log for details.'));
            else
                % Find most recent
                dates = {projFile.date};
                [time, mostRecent] = max(datenum(dates));
                projFile = fullfile(FileNames.filelocation,projFile(mostRecent).name);
                % Rename to correct naming convention
                movefile(projFile, fullfile(FileNames.filelocation, FileNames.eogProjFileName));
            end
        end
    else
        %         GUI.ErrorMessage(GUI.ErrorMessage.GENERIC_WARNING, sprintf('%s\n%s','Eye blink events not found.',' Cannot compute projections'));
    end
else
    disp('Eye blink SSP complete')
end
end

function nogui_remove_ecg_artifact(FileNames)
% remove_ecg_artifact: Get heartbeat events and create an SSP projection
%
% USAGE:    remove_ecg_artifact(mc, FileNames)
%
% INPUT:    mc = megclinic instance
%           FileNames = Filename structure
%
%
% Author: Elizabeth Bock, 2009
% --------------------------- Script History ------------------------------
% EB 27-AUG-2009  Creation
% -------------------------------------------------------------------------

% ------------- Read the raw file -----------------------------------------
% check to see if ecg event file already exists
checkfile = fullfile(FileNames.filelocation, FileNames.eventFileName);
fileExists = exist(checkfile,'file');
if (fileExists == 0)
    try
        [fiffsetup] = fiff_setup_read_raw(FileNames.filename);
    catch ME
        GUI.ErrorMessage(GUI.ErrorMessage.GENERIC_ERROR, ME.message);
        disp(ME.message)
        return
    end
    
    [ecg, channelType] = nogui_get_ecg(fiffsetup);
    
    if isempty(ecg)
        return;
    end
    
    Options.sampRate = fiffsetup.info.sfreq;
    firstSamp = fiffsetup.first_samp;
    
    % -------------- ECG -------------------------------------------------------
    % detect ecg events and write to .eve file
    
    if channelType %MEG
        Options.percentThresh = 90;     %(qrs detection threshold - percent)
        Options.noiseThresh = 2.5;           %(number of std from mean to include for detection)
        Options.maxCrossings = 3;         %(max number of crossings)
        
    else % ECG
        Options.percentThresh = 60;     %(qrs detection threshold - percent)
        Options.noiseThresh = 2.5;           %(number of std from mean to include for detection)
        Options.maxCrossings = 3;         %(max number of crossings)
    end
    
    Options.minBeats = 10;
    Options.ecgType = 999;
    
    %     mc.setMessage(GUI.Config.M_MAKE_ECG_EVENTS);
    ecg_events = nogui_qrsDet2(ecg, Options);
    % Use all ecg events for ssp and no longer limit the ssp to only 50 events. Used to be:    %% maxEvents = min(50, length(ecg_events));
    maxEvents = length(ecg_events);
    if ~isempty(ecg_events)
        nogui_writeEventFile(fullfile(FileNames.filelocation, FileNames.eventFileName), firstSamp, ecg_events(1:maxEvents), Options.ecgType);
        %         %         mc.setMessage(GUI.Config.M_ECG_EVENTS_WRITTEN);
    end
end

% -------------- SSP ------------------------------------------------------
% check for existing projection and event files
projExists = exist(fullfile(FileNames.filelocation, FileNames.ecgProjFileName),'file');
eveExists = exist(fullfile(FileNames.filelocation, FileNames.eventFileName), 'file');
if ~projExists
    if eveExists
        % create new ssp operators
        %         mc.setMessage(GUI.Config.M_MAKE_SSP);
        nogui_make_proj_operator(FileNames.filelocation, FileNames.filename, FileNames.eventFileName, 'ECG');
        
        projExists = exist(fullfile(FileNames.filelocation, FileNames.ecgProjFileName),'file');
        % If the proj file does not exist, then check to be sure MNE named
        % the file correctly.  If not, change the file name.
        if ~projExists
            % find all -proj files for ecg
            projFile = dir(fullfile(FileNames.filelocation, '*_ecg-proj.fif'));
            if isempty(projFile)
                % No proj files for ecg exist
                %                 GUI.ErrorMessage(GUI.ErrorMessage.GENERIC_WARNING, sprintf('%s\n%s','Heartbeat projections not created.',' See log for details'));
            else
                % Find most recent
                dates = {projFile.date};
                [time, mostRecent] = max(datenum(dates));
                projFile = fullfile(FileNames.filelocation,projFile(mostRecent).name);
                % Rename to correct naming convention
                movefile(projFile, fullfile(FileNames.filelocation, FileNames.ecgProjFileName));
            end
        end
    else
        %         GUI.ErrorMessage(GUI.ErrorMessage.GENERIC_WARNING, sprintf('%s\n%s','Heartbeat events not found.',' Cannot compute projections.'));
    end
    
else
    disp('Heartbeat SSP complete')
end

end

function [ecg, channelType] = nogui_get_ecg(fiffsetup)
% GET_ECG: Extract ECG from raw file
%
% USAGE:    ecg = get_ecg(fiffsetup);
%
% INPUT:    fiffsetup = fiff file info structure (mne)
%
% OUTPUT:   ecg = channel signal containing ecg information
%           channelType = 0, ECG channel
%           channelType = 1, MEG channel
%
% Author: Elizabeth Bock, 2009
% --------------------------- Script History ------------------------------
% EB 27-AUG-2009  Creation
% -------------------------------------------------------------------------

% config = upper(char(ArtifactClean.CleanConfig.ECG_CHAN));

config = 'AUTO';

% ----------------------- Determine Channel for ECG information -----------
if strcmp(config, 'AUTO') || strcmp(config, '')
    channel = 'ECG';
    channelType = 0; % this is default for ECG channel
else
    channel = config;
    channelType = 1;
end


channelNames = fiffsetup.info.ch_names;
ch_ECG = strmatch(channel,channelNames);

if isempty(ch_ECG) % ECG channel does not exist
    [ch_ECG,okButton] = listdlg('Name', 'No ECG Channel', 'PromptString', 'Select Another Channel?:',...
        'SelectionMode','single',...
        'ListString',channelNames); % prompt user to choose an MEG channel
    if okButton
        channelType = 1; % MEG channel choosen for analysis
    else
        GUI.ErrorMessage(GUI.ErrorMessage.GENERIC_ERROR, 'No ECG Channel Available for Cleaning.');
        ecg = [];
        return;
    end
end

% ------------------ Get channel signal from raw file ---------------------
start_samp = fiffsetup.first_samp;
end_samp = fiffsetup.last_samp;
[ecg,ecgtimes] = fiff_read_raw_segment(fiffsetup, start_samp ,end_samp, ch_ECG);

end

function clean_events = nogui_qrsDet2(ecg, Options)
% qrsDet2: Detect QRS events from ECG signal
%
% USAGE: clean_events = qrsDet2(ecg, Options)
%
% INPUT:    ecg = channel signal containing ECG information
%
%           Options.sampRate        sampling rate of the signal
%           Options.percentThresh   qrs detection threshold (percent)
%           Options.noiseThresh     number of std from mean to include for detection)
%           Options.maxCrossings    max number of crossings
%           Options.minBeats        minimum number of heartbeats needed for event file
%
% OUTPUT:   clean_events = qrs event times
%
% Author: Elizabeth Bock, 2009
% --------------------------- File History ------------------------------
% EB 27-AUG-2009  Creation
% -------------------------------------------------------------------------

% logFile = GUI.MCLogFile;
clean_events = [];

% --------------------- Set up blanking period ----------------------------
% No events can be detected during the blanking
% period.  This blanking period assumes no heartrate will be faster than
% 120 bpm.
minInterval = round((60*Options.sampRate)/120);
BLANK_PERIOD = minInterval;

% ------------------ Determine the ecg channel ----------------------------
ecgChan = 1;
[row, col] = size(ecg);
% If there is more than one channel labeled as ECG, let the user choose.
if row > 1
    %     logFile.write('2 ECG Leads Detected')
    % plot 10 secs of each lead
    subplot(2,1,1)
    plot(ecg(1,1:10000))
    hold on
    subplot(2,1,2)
    plot(ecg(2,1:10000))
    
    disp('Two ECG leads were detected, which would you like to use for analysis?')
    reply = input('1 or 2?');
    if reply == 2
        ecgChan = 2;
    else
        ecgChan = 1;
    end
end

% -------------------------- Bandpass filter the signal--------------------
filtecg = bst_bandpass_fft(ecg(ecgChan,:), Options.sampRate, 5, 35);
lenpts = length(filtecg);
% Absolute value
absecg = abs(filtecg);

% -------------------------- Determine threshold --------------------------
init = round(Options.sampRate); % One second
% Find the average of the maximum of each of the first three seconds
maxpt(1) = max(absecg(1:init));
maxpt(2) = max(absecg(init:2*init));
maxpt(3) = max(absecg(init*2:init*3));
init_max = mean(maxpt);
% Threshold is percent of max
thresh_value = Options.percentThresh/100;
qrs_event.thresh1 = init_max*thresh_value;
if qrs_event.thresh1 < 0.02e-3
    return;
end

% ------------------------- Find events -----------------------------------
qrs_event.filtecg = filtecg;
qrs_event.time=[];
k=1;
i=1;

while i < lenpts-BLANK_PERIOD+1
    if absecg(i) > qrs_event.thresh1 % signal exceeds thresh
        window = absecg(i:i+BLANK_PERIOD);
        [maxPoint, maxTime] = max(window(1:BLANK_PERIOD/2)); % max of the window
        rms = sqrt(mean(window.^2)); % rms of signal for this window
        
        % Find the number of threshold crossings in the window
        x=find(window>qrs_event.thresh1);
        y=diff(x);
        numcross = length(find(y>1))+1;
        
        qrs_event.time(k) = maxTime+i; % time of max value
        qrs_event.ampl(k) = maxPoint; % max value
        qrs_event.numcross(k) = numcross;
        qrs_event.rms(k) = rms;
        
        i=i+BLANK_PERIOD; % skip ahead past blank period
        k=k+1; % increment event count
    else
        i=i+1; % increment to next point
    end
end

if ~isempty(qrs_event.time)
    % Exclude events that do not meet noise criteria
    rms_mean = mean(qrs_event.rms); % mean rms of all events
    rms_std = std(qrs_event.rms); % std of rms for all events
    rms_thresh = rms_mean+(rms_std*Options.noiseThresh); % rms threshold
    b = find(qrs_event.rms < rms_thresh); % find events less than rms threshold
    a = qrs_event.numcross(b);
    c = a < Options.maxCrossings; % find events with threshold crosses less than desired value
    clean_events = qrs_event.time(b(c));
end

nBeats = length(clean_events);
% logFile.write(['Only ' num2str(nBeats) ' heart beats were detected.  This recording will not be cleaned.']);
if nBeats < Options.minBeats
    GUI.ErrorMessage(GUI.ErrorMessage.GENERIC_WARNING, ['Only ' num2str(nBeats) ' heart beats were detected.  This recording will not be cleaned.']);
    clean_events = [];
end
end

function nogui_writeEventFile(eventFileName, firstSamp, events, eventType)
% WRITEEVENTFILE: creates .eve text file (MNE) from a list of events
%
% USAGE:  writeEventFile(fid, firstSamp, lastSamp, sampRate, events, eventType);
%
% INPUT:
%   - eventFileName : name of file to save event
%   - firstSamp   : first valid sample of raw data (type double)
%   - event       : array of events (samples)
%   - eventType   : type of events (ie 1000)
% OUTPUT:
%   _sss-eve.fif file with events (MNE format)
%
% Author: Elizabeth Bock, 2009
% --------------------------- Script History ------------------------------
% EB 27-AUG-2009  Creation
% -------------------------------------------------------------------------

% logFile = GUI.MCLogFile;
% Adjust the event samples for MNE
firstSamp = double(firstSamp);
adj_events = events+firstSamp;

% Include psuedo event
eventlist(1,1) = firstSamp;
eventlist(1,2) = 0;
eventlist(1,3) = 0;

% Include all other events
num_events = length(events);
% logFile.write(['ECG Events Detected:' char(num_events)]);
eventlist(2:num_events+1,1) = adj_events;
eventlist(2:num_events+1,2) = zeros(1,num_events);
eventlist(2:num_events+1,3) = ones(1,num_events)*eventType;

% Write to .fif file
mne_write_events(eventFileName,eventlist);

end

function make_megclinic_xml_templates(subjIDpath, outDirpath)
% Create session.xml and raw.xml in the run directory
subjIDpath = char(subjIDpath);
outDirpath = char(outDirpath);

% In your example, run folder is subjIDpath (e.g., .../260202)
runDir = subjIDpath;

analysesDir = fullfile(runDir, 'analyses');
reportDir   = fullfile(runDir, 'report');
if exist(analysesDir,'dir') ~= 7, mkdir(analysesDir); end
if exist(reportDir,'dir') ~= 7, mkdir(reportDir); end

sessionFile = fullfile(runDir, 'session.xml');
rawFile     = fullfile(runDir, 'raw.xml');

% ---- session.xml
fid = fopen(sessionFile, 'w');
if fid < 0, error('Cannot write: %s', sessionFile); end
fprintf(fid, '<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid, '<SESSION><TYPE>');
fprintf(fid, '<RAWDIR>%s</RAWDIR>', runDir);
fprintf(fid, '<MRIDIR/><EMPTY/><INTERICTALSPIKES/><INTERICTALFILES/><SPONT/>');
fprintf(fid, '<FUNCTIONALEVENTS/><FUNCTIONALFILES/><PROBEPOSITION/>');
fprintf(fid, '<BSTDIR/><BSTPROTOCOL/><BSTPROTOCOLTYPE/><BSTMRI/><BSTPIPELINE/>');
fprintf(fid, '<ANALYSIS>%s</ANALYSIS>', analysesDir);
fprintf(fid, '<REPORT>%s</REPORT>', reportDir);
fprintf(fid, '<FREQBANDS>delta: 1, 4; theta: 4, 8; alpha: 8, 12; beta: 13, 35; gamma1: 40, 80; gamma2: 80, 150; gamma3: 150, 250; gamma4: 250, 350; gamma5: 350, 550;</FREQBANDS>');
fprintf(fid, '<FXNLIMPORTSTRT/><FXNLIMPORTSTOP/><FXNLBASESTRT/><FXNLBASESTOP/>');
fprintf(fid, '</TYPE></SESSION>\n');
fclose(fid);

% ---- raw.xml
fid = fopen(rawFile, 'w');
if fid < 0, error('Cannot write: %s', rawFile); end
fprintf(fid, '<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid, '<WORKFLOW><TYPE>');
fprintf(fid, '<RUNDIR>%s</RUNDIR>', runDir);
fprintf(fid, '<RAW>raw</RAW>');
fprintf(fid, '<AVE/><RAWSSS/><AVESSS/><ECGEVE/><ECGPROJ/><ECGCLEAN/>');
fprintf(fid, '<ARTEVE/><ARTPROJ/><ARTCLEAN/><EOGEVE/><EOGPROJ/><EOGCLEAN/>');
fprintf(fid, '<ONGNGPROJ/><ONGNGCLEAN/><EVENTLIST/><FXNLEVENTS/><AVEDESC/><AVECLEAN/>');
fprintf(fid, '<MASK/><STIMSOURCE/><USEREVE/><AVEUSEREVE/><CUSTOMAVEDESC/>');
fprintf(fid, '<MRIDIR/><MRIFIFF/><XFITDIR/><CFIT/><DIPOLES/><DIPOLEIMG/>');
fprintf(fid, '<BSTDIR/><BSTPROTOCOL/><BSTPROTOCOLTYPE/><BSTMRI/><BSTMEG/>');
fprintf(fid, '<REPORT/><XML/><EMTPY/><PROBEPOSITION/>');
fprintf(fid, '<ANALYSES>%s</ANALYSES>', analysesDir);
fprintf(fid, '</TYPE></WORKFLOW>\n');
fclose(fid);

fprintf('Wrote XML templates:\n  %s\n  %s\n', sessionFile, rawFile);
end
