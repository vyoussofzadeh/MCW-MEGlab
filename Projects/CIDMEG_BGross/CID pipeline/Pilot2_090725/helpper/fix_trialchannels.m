
cd('/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/Brainstorm_db/CID/data/mcwa086_v1/mcwa086_v1_Run_03')


% fix_bst_trial_channel_mismatch.m
% Run this in the folder containing your Brainstorm trial files (data_*_trial*.mat)

% fix_bst_trial_channel_mismatch.m
% Run this in the folder containing your Brainstorm trial files (data_*_trial*.mat)

% clear; 
clc;

files = dir(fullfile(pwd, 'data_*_trial*.mat'));
timestamp = datestr(now,'yyyymmdd_HHMMSS');
backupDir = fullfile(pwd, ['backup_before_fix_' timestamp]);
if ~exist(backupDir, 'dir'), mkdir(backupDir); end

nPatched = 0; nOK = 0; nSkipped = 0;

for k = 1:numel(files)
    fn = fullfile(files(k).folder, files(k).name);
    s = load(fn);  % load all vars (F, ChannelFlag, Time, etc.)
    
%     if ~isfield(s,'F') || ~isfield(s,'ChannelFlag')
%         fprintf('SKIP (no F/ChannelFlag): %s\n', files(k).name);
%         nSkipped = nSkipped + 1;
%         continue;
%     end
    
    % Ensure column vector
    s.ChannelFlag = s.ChannelFlag(:);
    
    nF = size(s.F, 1);
    nFlag = numel(s.ChannelFlag);
    
    if nF == nFlag
        fprintf('OK   (%3d ch): %s\n', nF, files(k).name);
        nOK = nOK + 1;
        continue;
    end
    
    % Backup original
    copyfile(fn, fullfile(backupDir, files(k).name));
    
    action = '';
    if nF < nFlag
        % Common case after adding 2 MISC channels to the channel file:
        % pad F with zeros at the end to match ChannelFlag length
        d = nFlag - nF;
        s.F = [s.F; zeros(d, size(s.F,2), class(s.F))];
        action = sprintf('Padded F with %d zero rows to match ChannelFlag (%d).', d, nFlag);
        
    else % nF > nFlag
        % Less common: more rows in F than flags. Extend ChannelFlag with "good" (1).
        d = nF - nFlag;
        s.ChannelFlag = [s.ChannelFlag; ones(d,1)];
        action = sprintf('Extended ChannelFlag by %d entries (set to 1) to match F (%d).', d, nF);
    end
    
    % Append to History if present
    if isfield(s,'History') && iscell(s.History)
        try
            s.History(end+1,1:3) = {datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
                'fix_bst_trial_channel_mismatch', action};
        catch
            s.History{end+1,1} = action;  % fallback
        end
    end
    
    save(fn, '-struct', 's');  % overwrite with fixed struct
    fprintf('FIX  %s | %s\n', files(k).name, action);
    nPatched = nPatched + 1;
end

fprintf('\nPatched: %d | OK: %d | Skipped: %d\nBackups: %s\n', ...
    nPatched, nOK, nSkipped, backupDir);
