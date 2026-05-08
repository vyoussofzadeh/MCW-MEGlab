dbRoot = '/MEG_data/epilepsy/*/brainstorm_db/data/';

files = dir(fullfile(dbRoot, 'protocol.mat'));

% Output containers
ProtocolFolder = {};
ProtocolFile   = {};
SubjectFolder  = {};
DbVersion      = {};
Status         = {};
ErrorMessage   = {};

for i = 1:numel(files)
    f = fullfile(files(i).folder, files(i).name);

    try
        vars = who('-file', f);

        % Try to extract subject/project folder name after /epilepsy/
        parts = split(files(i).folder, filesep);
        idx = find(strcmp(parts, 'epilepsy'), 1, 'last');

        if ~isempty(idx) && numel(parts) >= idx + 1
            subj = parts{idx + 1};
        else
            subj = '';
        end

        if ismember('DbVersion', vars)
            S = load(f, 'DbVersion');

            ProtocolFolder{end+1,1} = files(i).folder;
            ProtocolFile{end+1,1}   = f;
            SubjectFolder{end+1,1}  = subj;
            DbVersion{end+1,1}      = string(S.DbVersion);
            Status{end+1,1}         = 'DbVersion found';
            ErrorMessage{end+1,1}   = '';

            fprintf('\nDbVersion found:\n');
            fprintf('Subject/Folder: %s\n', subj);
            fprintf('File: %s\n', f);
            disp(S.DbVersion)

        else
            ProtocolFolder{end+1,1} = files(i).folder;
            ProtocolFile{end+1,1}   = f;
            SubjectFolder{end+1,1}  = subj;
            DbVersion{end+1,1}      = '';
            Status{end+1,1}         = 'DbVersion not found';
            ErrorMessage{end+1,1}   = '';
        end

    catch ME
        ProtocolFolder{end+1,1} = files(i).folder;
        ProtocolFile{end+1,1}   = f;
        SubjectFolder{end+1,1}  = '';
        DbVersion{end+1,1}      = '';
        Status{end+1,1}         = 'Could not read file';
        ErrorMessage{end+1,1}   = ME.message;

        fprintf('Could not read: %s\nReason: %s\n', f, ME.message);
    end
end

% Create table
T = table(SubjectFolder, ProtocolFolder, ProtocolFile, DbVersion, Status, ErrorMessage);

% Display table
disp(T)

% Save table
outCsv = '/MEG_data/epilepsy/Brainstorm_DbVersion_report.csv';
outMat = '/MEG_data/epilepsy/Brainstorm_DbVersion_report.mat';

% writetable(T, outCsv);
% save(outMat, 'T');

fprintf('\nSaved table to:\n%s\n%s\n', outCsv, outMat);

%%

dbRoot = '/MEG_data/epilepsy/*/brainstorm_db/data/';
refVersion = 5.0300;   % BS2020-compatible reference version

files = dir(fullfile(dbRoot, 'protocol.mat'));

SubjectFolder  = {};
ProtocolFolder = {};
ProtocolFile   = {};
DbVersion      = [];

for i = 1:numel(files)

    f = fullfile(files(i).folder, files(i).name);

    try
        vars = who('-file', f);

        if ismember('DbVersion', vars)
            S = load(f, 'DbVersion');

            % Convert DbVersion to numeric if needed
            if isnumeric(S.DbVersion)
                thisVersion = double(S.DbVersion);
            else
                thisVersion = str2double(string(S.DbVersion));
            end

            % Keep only databases with higher DbVersion
            if thisVersion == refVersion

                % Extract epilepsy subject/case folder
                parts = split(files(i).folder, filesep);
                idx = find(strcmp(parts, 'epilepsy'), 1, 'last');

                if ~isempty(idx) && numel(parts) >= idx + 1
                    subj = parts{idx + 1};
                else
                    subj = '';
                end

                SubjectFolder{end+1,1}  = subj;
                ProtocolFolder{end+1,1} = files(i).folder;
                ProtocolFile{end+1,1}   = f;
                DbVersion(end+1,1)      = thisVersion;
            end
        end

    catch ME
        fprintf('Could not read: %s\nReason: %s\n', f, ME.message);
    end
end

% Create and sort table
T_high = table(SubjectFolder, ProtocolFolder, ProtocolFile, DbVersion);
T_high = sortrows(T_high, 'DbVersion', 'descend');

% Display result
disp(T_high)

% % Save report
% outCsv = '/MEG_data/epilepsy/Brainstorm_high_DbVersion_report.csv';
% outMat = '/MEG_data/epilepsy/Brainstorm_high_DbVersion_report.mat';
% 
% writetable(T_high, outCsv);
% save(outMat, 'T_high');
% 
% fprintf('\nFound %d database(s) with DbVersion > %.4f\n', height(T_high), refVersion);
% fprintf('Saved report to:\n%s\n%s\n', outCsv, outMat);

%%
% Make sure Brainstorm is closed before running this

newDbVersion = 5.0100;
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

for i = 1:height(T_high)

    f = T_high.ProtocolFile{i};

    fprintf('\nProcessing:\n%s\n', f);

    try
        % Load current DbVersion
        S = load(f, 'DbVersion');

        fprintf('Current DbVersion: %.4f\n', S.DbVersion);

        % Create backup copy first
        backupFile = fullfile(fileparts(f), ...
            ['protocol_backup_DbVersion_' num2str(S.DbVersion, '%.4f') '_' timestamp '.mat']);

        copyfile(f, backupFile);
        fprintf('Backup saved:\n%s\n', backupFile);

        % Update DbVersion
        DbVersion = newDbVersion;
        save(f, 'DbVersion', '-append');

        % Confirm update
        Scheck = load(f, 'DbVersion');
        fprintf('Updated DbVersion: %.4f\n', Scheck.DbVersion);

    catch ME
        fprintf('Could not update:\n%s\nReason: %s\n', f, ME.message);
    end
end