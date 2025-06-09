function ft_raw = do_read_neuromag_raw_nofilter(datafile)
% READ_NEUROMAG_RAW_NOFILTER  Load a .fif that still has MaxShield on it
%                             and return a FieldTrip raw struct.
%
%  ft_raw = read_neuromag_raw_nofilter('/path/to/file_raw.fif')

% -------------------------------------------------------------------------
% 0)  Make sure Brainstorm is running and on the path
% -------------------------------------------------------------------------
if ~exist('brainstorm', 'file')
    addpath('/usr/local/MATLAB_Tools/BS_2024');   % adjust to your install
end
% if ~brainstorm('status')
brainstorm nogui
% brainstorm
% end

% -------------------------------------------------------------------------
% 1)  Re-use (or create) a throw-away protocol
% -------------------------------------------------------------------------
protoName = 'TMP_MaxShield';
iProto    = bst_get('Protocol', protoName);
if isempty(iProto)
    gui_brainstorm('CreateProtocol', protoName, 0, 0);  % default anatomy
else
    gui_brainstorm('SetCurrentProtocol', iProto);
end

% -------------------------------------------------------------------------
% 2)  Re-use (or create) a temp subject
% -------------------------------------------------------------------------
subjName = 'TmpSubj';
if isempty(bst_get('Subject', subjName))
    db_add_subject(subjName);
end

%%

disp('DO THE MANUAL IMPORT, if needed')
disp('hit enter')
% pause,

%%

% -------------------------------------------------------------------------
% 3)  Import the raw file as a *link* (no SSS/tSSS)
% -------------------------------------------------------------------------
sFiles = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
    'subjectname',    'TmpSubj', ...
    'datafile',       {datafile, 'FIF'}, ...
    'channelreplace', 0, ...
    'channelalign',   1, ...
    'evtmode',        'value');

% -------------------------------------------------------------------------
% 4)  Convert the link to FieldTrip format
% -------------------------------------------------------------------------
protInfo = bst_get('ProtocolInfo');
dataRoot = protInfo.STUDIES;

linkFile = fullfile(dataRoot, sFiles.FileName);
chanFile = fullfile(dataRoot, sFiles.ChannelFile);

[ft_raw, D] = out_fieldtrip_data(linkFile, chanFile, 'MEG', 0);

% add FieldTrip bookkeeping fields
ft_raw.time       = {D.Time};
ft_raw.fsample    = 1 / mean(diff(D.Time));
ft_raw.sampleinfo = [1 numel(D.Time)];
ft_raw.grad.chantype = ft_raw.grad.chantype(:);   % column, not row

% Optional: tidy the temporary files to keep the DB small
% bst_delete.File(sFiles);
end
