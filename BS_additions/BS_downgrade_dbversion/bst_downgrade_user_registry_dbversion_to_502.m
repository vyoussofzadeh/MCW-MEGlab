% Local Brainstorm registry file
bstUserFile = fullfile(getenv('HOME'), '.brainstorm', 'brainstorm.mat');

% Target BS2020-compatible version
newDbVersion = 5.0200;

% Check file exists
if ~exist(bstUserFile, 'file')
    error('Could not find: %s', bstUserFile);
end

% Load current DbVersion
S = load(bstUserFile, 'DbVersion');

fprintf('Current DbVersion: %.4f\n', S.DbVersion);

% Make backup first
backupFile = [bstUserFile '.backup_' datestr(now, 'yyyymmdd_HHMMSS')];
copyfile(bstUserFile, backupFile);

fprintf('Backup saved:\n%s\n', backupFile);

% Update DbVersion only
DbVersion = newDbVersion;
save(bstUserFile, 'DbVersion', '-append');

% Confirm
Scheck = load(bstUserFile, 'DbVersion');
fprintf('Updated DbVersion: %.4f\n', Scheck.DbVersion);