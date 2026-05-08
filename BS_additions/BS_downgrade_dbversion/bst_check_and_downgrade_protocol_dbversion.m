bstUserFile = '/home/vyoussof/.brainstorm/brainstorm.mat';

a = load(bstUserFile);
P = a.ProtocolsListInfo;

ProtocolName   = {};
SubjectDir     = {};
StudyDir       = {};
ProtocolMat    = {};
DbVersion      = [];
Status         = {};
ErrorMessage   = {};

for i = 1:numel(P)

    % Protocol name/comment
    if isfield(P, 'Comment')
        protocolName = P(i).Comment;
    else
        protocolName = ['Protocol_' num2str(i)];
    end

    % Try to get STUDIES and SUBJECTS folders
    studyDir = '';
    subjDir  = '';

    if isfield(P, 'STUDIES')
        studyDir = P(i).STUDIES;
    end

    if isfield(P, 'SUBJECTS')
        subjDir = P(i).SUBJECTS;
    end

    % Brainstorm protocol.mat is usually in the STUDIES/data folder
    candidateFiles = {};

    if ~isempty(studyDir)
        candidateFiles{end+1} = fullfile(studyDir, 'protocol.mat');
    end

    if ~isempty(subjDir)
        candidateFiles{end+1} = fullfile(fileparts(subjDir), 'data', 'protocol.mat');
    end

    protocolFile = '';
    for j = 1:numel(candidateFiles)
        if exist(candidateFiles{j}, 'file')
            protocolFile = candidateFiles{j};
            break
        end
    end

    try
        if isempty(protocolFile)
            thisVersion = NaN;
            statusMsg = 'protocol.mat not found';
            errMsg = '';
        else
            vars = who('-file', protocolFile);

            if ismember('DbVersion', vars)
                S = load(protocolFile, 'DbVersion');
                thisVersion = double(S.DbVersion);
                statusMsg = 'DbVersion found';
                errMsg = '';
            else
                thisVersion = NaN;
                statusMsg = 'DbVersion not found';
                errMsg = '';
            end
        end

    catch ME
        thisVersion = NaN;
        statusMsg = 'Could not read protocol.mat';
        errMsg = ME.message;
    end

    ProtocolName{end+1,1} = protocolName;
    SubjectDir{end+1,1}   = subjDir;
    StudyDir{end+1,1}     = studyDir;
    ProtocolMat{end+1,1}  = protocolFile;
    DbVersion(end+1,1)    = thisVersion;
    Status{end+1,1}       = statusMsg;
    ErrorMessage{end+1,1} = errMsg;
end

T_protocols = table(ProtocolName, SubjectDir, StudyDir, ProtocolMat, DbVersion, Status, ErrorMessage);

disp(T_protocols)

%%
T_high = T_protocols(T_protocols.DbVersion > 5.0200, :);
disp(T_high)

%%
% Make sure Brainstorm is closed before running this

newDbVersion = 5.0200;
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

for i = 1:height(T_high)

    f = T_high.ProtocolMat{i};

    fprintf('\nProcessing:\n%s\n', f);

    try
        % Load current DbVersion
        S = load(f, 'DbVersion');
        fprintf('Current DbVersion: %.4f\n', S.DbVersion);

        % Make backup copy
        backupFile = fullfile(fileparts(f), ...
            ['protocol_backup_DbVersion_' num2str(S.DbVersion, '%.4f') '_' timestamp '.mat']);

        copyfile(f, backupFile);
        fprintf('Backup saved:\n%s\n', backupFile);

        % Update DbVersion only
        DbVersion = newDbVersion;
        save(f, 'DbVersion', '-append');

        % Confirm
        Scheck = load(f, 'DbVersion');
        fprintf('Updated DbVersion: %.4f\n', Scheck.DbVersion);

    catch ME
        fprintf('Could not update:\n%s\nReason: %s\n', f, ME.message);
    end
end