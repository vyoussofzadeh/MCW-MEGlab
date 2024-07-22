function ecp_import_raw(RawFiles, iSubject)

%%
ImportOptions = [];
ImportOptions.UseEvents = 0;
ImportOptions.TimeRange = [];
ImportOptions.EventsTimeRange= [-0.1000    0.3000];
ImportOptions.GetAllEpochs= 0;
ImportOptions.iEpochs= 1;
ImportOptions.SplitRaw= 0;
ImportOptions.SplitLength= 2;
ImportOptions.Resample= 0;
ImportOptions.ResampleFreq= 0;
ImportOptions.UseCtfComp= 1;
ImportOptions.UseSsp= 1;
ImportOptions.RemoveBaseline= 'no';
ImportOptions.BaselineRange= [];
ImportOptions.BaselineSensorType= '';
ImportOptions.events= [];
ImportOptions.CreateConditions= 0;
ImportOptions.ChannelReplace= 1;
ImportOptions.ChannelAlign= 1;
ImportOptions.IgnoreShortEpochs= 1;
ImportOptions.EventsMode= 'ask';
ImportOptions.EventsTrackMode= 'ask';
ImportOptions.EventsTypes= '';
ImportOptions.DisplayMessages= 1;
ImportOptions.Precision= [];

[sFile, ChannelMat] = in_fopen_fif_edit(RawFiles, ImportOptions);

sFile.events = struct_fix_events(sFile.events);

% === SORT BY NAME ===
% Remove the common components
[tmp__, evtLabels] = str_common_path({sFile.events.label});
% Try to convert all the names to numbers
evtNumber = cellfun(@str2num, evtLabels, 'UniformOutput', 0);

% === ADD COLOR ===
if isempty(sFile.events(1).color)
    ColorTable = panel_record('GetEventColorTable');
    for i = 1:length(sFile.events)
        iColor = mod(i-1, length(ColorTable)) + 1;
        sFile.events(i).color = ColorTable(iColor,:);
    end
end

[fPath, fBase] = bst_fileparts(RawFiles);

nSes = 1; iFile = 1; FileFormat = 'FIF';
[fPath, fBase] = bst_fileparts(RawFiles);
sSubject = bst_get('Subject', iSubject, 1);
ChannelMat = bst_history('add', ChannelMat, 'import', ['Link to file: ' RawFiles ' (Format: ' FileFormat ')']);
% Remove fiducials only from polhemus and ascii files
isRemoveFid = ismember(FileFormat, {'MEGDRAW', 'POLHEMUS', 'ASCII_XYZ', 'ASCII_NXYZ', 'ASCII_XYZN', 'ASCII_XYZ_MNI', 'ASCII_NXYZ_MNI', 'ASCII_XYZN_MNI', 'ASCII_NXY', 'ASCII_XY', 'ASCII_NTP', 'ASCII_TP'});
% Perform the NAS/LPA/RPA registration for some specific file formats
isAlign = ismember(FileFormat, {'NIRS-BRS','NIRS-SNIRF'});
% Detect auxiliary EEG channels
ChannelMat = channel_detect_type(ChannelMat, isAlign, isRemoveFid);
% Do not align data coming from Brainstorm exported files (already aligned)
ConditionName = ['@raw' fBase];
NewComment = 'Link to raw file';
[sExistStudy, iExistStudy] = bst_get('StudyWithCondition', bst_fullfile(sSubject.Name, file_standardize(ConditionName, 1)));
studyDate = sFile.acq_date;
iOutputStudy = db_add_condition(sSubject.Name, ConditionName, [], studyDate);
% Get output study
sOutputStudy = bst_get('Study', iOutputStudy);
% Get the study in which the channel file has to be saved
[sChannel, iChannelStudy] = bst_get('ChannelForStudy', iOutputStudy);

Tolerance = [];
[ChannelFile, ChannelMat, ImportOptions.ChannelReplace, ImportOptions.ChannelAlign, Modality, Tolerance] = db_set_channel_edit(iChannelStudy, ChannelMat, ImportOptions.ChannelReplace, ImportOptions.ChannelAlign, Tolerance);
% [ChannelFile, ChannelMat, ImportOptions.ChannelReplace, ImportOptions.ChannelAlign, Modality, Tolerance] = db_set_channel(iChannelStudy, ChannelMat, ImportOptions.ChannelReplace, ImportOptions.ChannelAlign, Tolerance);

sFileOut = sFile;


ProtocolInfo = bst_get('ProtocolInfo');


% ===== SAVE LINK FILE =====
% Build output filename
NewBstFile = bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(sOutputStudy.FileName), ['data_0raw_' fBase '.mat']);
% Build output structure
NewMat = db_template('DataMat');
NewMat.F           = sFileOut;
NewMat.Comment     = NewComment;
NewMat.ChannelFlag = sFileOut.channelflag;
NewMat.Time        = sFileOut.prop.times;
NewMat.DataType    = 'raw';
NewMat.Device      = sFileOut.device;
% Compumedics: add start time to the file comment
if strcmpi(sFileOut.format, 'EEG-COMPUMEDICS-PFS') && isfield(sFileOut.header, 'rda_startstr') && ~isempty(sFileOut.header.rda_startstr)
    NewMat.Comment = [NewMat.Comment ' [' sFileOut.header.rda_startstr ']'];
end
% Add history field
NewMat = bst_history('add', NewMat, 'import', ['Link to raw file: ' RawFiles]);
% Save file on hard drive
bst_save(NewBstFile, NewMat, 'v6');
% Add file to database
sOutputStudy = db_add_data(iOutputStudy, NewBstFile, NewMat);
% Return new file
% OutputFiles{end+1} = NewBstFile;
% ===== UPDATE DATABASE =====
% Update links
db_links('Study', iOutputStudy);
% Refresh both data node and channel node
iUpdateStudies = unique([iOutputStudy, iChannelStudy]);
panel_protocols('UpdateNode', 'Study', iUpdateStudies);
% ===== ADD SYNCHRONIZED VIDEOS =====
% Look for video files with the same name, add them to the raw folder if found
for fExt = {'.avi','.AVI','.mpg','.MPG','.mpeg','.MPEG','.mp4','.MP4','.mp2','.MP2','.mkv','.MKV','.wmv','.WMV','.divx','.DIVX','.mov','.MOV'}
    VideoFile = bst_fullfile(fPath, [fBase, fExt{1}]);
    if file_exist(VideoFile)
        import_video(iOutputStudy, VideoFile);
        break;
    end
end
% If something was imported
if ~isempty(iOutputStudy)
    % Select the data study node
    panel_protocols('SelectStudyNode', iOutputStudy);
    % Save database
    db_save();
end
%         if isProgress
bst_progress('stop');
%         end
end