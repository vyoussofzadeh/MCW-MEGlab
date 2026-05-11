function trials = build_trials_from_mat(filesS, filesN, varargin)
% Build aligned MEG epoch cells + labels + groups from annotated .mat files.
% trials.data   : {N×1} cell of [C×T] epochs
% trials.labels : N×1 categorical ('NoSpike','Spike')
% trials.groups : {N×1} group ids (subject/session)
% trials.refLbl : C×1 label list
%
% Optional name/value:
%   'Verbose'   (true/false) default true
%   'UseWaitbar'(true/false) default false  (forces MATLAB waitbar even if ft_progress exists)

p = inputParser;
addParameter(p,'Verbose',true,@islogical);
addParameter(p,'UseWaitbar',false,@islogical);
parse(p,varargin{:});
VERBOSE    = p.Results.Verbose;
USEWAITBAR = p.Results.UseWaitbar;

nS = numel(filesS); nN = numel(filesN);
if VERBOSE
    fprintf('[build_trials] Spike files: %d | NoSpike files: %d\n', nS, nN);
end

% ---- 1) Determine reference label order from the first Spike file ----
refLbl = [];
for k = 1:nS
    S = load(filesS(k).name);
    if isfield(S,'anot_data_all') && ~isempty(S.anot_data_all)
        D0 = S.anot_data_all{1};
        allLbl  = D0.label(:);
        megMask = startsWith(allLbl,'MEG');
        refLbl  = allLbl(megMask);
        if VERBOSE
            fprintf('[build_trials] Reference labels taken from %s (C=%d)\n', filesS(k).name, numel(refLbl));
        end
        break
    end
    if VERBOSE, fprintf('[build_trials] %s has no anot_data_all, skipping.\n', filesS(k).name); end
end
assert(~isempty(refLbl),'Could not determine reference MEG labels from spike files.');

% ---- 2) Read files into epoch cells (Spike then NoSpike) with progress ----
[dataS, groupsS] = read_filelist(filesS, refLbl, 'Label','Spike',   'Verbose',VERBOSE,'UseWaitbar',USEWAITBAR);
[dataN, groupsN] = read_filelist(filesN, refLbl, 'Label','NoSpike', 'Verbose',VERBOSE,'UseWaitbar',USEWAITBAR);

% ---- 3) Remove any NaN epochs class-wise (rare but safer) ----
isBad = @(x) any(isnan(x(:)));
keepS = ~cellfun(isBad, dataS);
keepN = ~cellfun(isBad, dataN);
if VERBOSE
    fprintf('[build_trials] Removed NaN epochs: Spike %d ? %d | NoSpike %d ? %d\n', ...
        numel(dataS), sum(keepS), numel(dataN), sum(keepN));
end
dataS = dataS(keepS); groupsS = groupsS(keepS);
dataN = dataN(keepN); groupsN = groupsN(keepN);

% ---- 4) Concatenate with categorical labels (fixed order) ----
data   = [dataN; dataS];
labels = categorical([zeros(numel(dataN),1); ones(numel(dataS),1)], [0 1], {'NoSpike','Spike'});
groups = [groupsN; groupsS];

% ---- 5) Sanity checks ----
chCounts = cellfun(@(x) size(x,1), data);
assert(all(chCounts==numel(refLbl)), 'Channel mismatch: not all epochs match refLbl length.');

% ---- 6) Return ----
trials = struct('data',{data}, 'labels',labels, 'groups',{groups}, 'refLbl',{refLbl});

if VERBOSE
    fprintf('[build_trials] DONE. Total epochs: %d  (NoSpike=%d, Spike=%d)\n', numel(data), numel(dataN), numel(dataS));
end
end


function [cells, groups] = read_filelist(filelist, refLbl, varargin)
% Load all anot_data_all epochs from a list of files, align to refLbl, with progress.
p = inputParser;
addParameter(p,'Label','',@ischar);
addParameter(p,'Verbose',true,@islogical);
addParameter(p,'UseWaitbar',false,@islogical);
parse(p,varargin{:});
tag       = p.Results.Label;
VERBOSE   = p.Results.Verbose;
USEWAIT   = p.Results.UseWaitbar;

cells  = {};
groups = {};
nFiles = numel(filelist);
totEpochs = 0;

% Choose progress mechanism
useFT   = ~USEWAIT && exist('ft_progress','file')==2;
if useFT
    ft_progress('init','text', '[%s] Reading %d files...', tag, nFiles);
    c = onCleanup(@() ft_progress('close'));
elseif USEWAIT
    wb = waitbar(0, sprintf('[%s] Reading 0/%d files...', tag, nFiles));
    c = onCleanup(@() (ishandle(wb) && close(wb)));
else
    c = onCleanup(@() 1);
end

for k = 1:nFiles
    S = load(filelist(k).name);
    if ~isfield(S,'anot_data_all') || isempty(S.anot_data_all)
        if VERBOSE, fprintf('[%s] %s: no anot_data_all, skipped.\n', tag, filelist(k).name); end
        prog(k,nFiles,totEpochs); %#ok<NASGU>
        continue;
    end

    gid = local_group_id(filelist(k).name);
    nThis = 0;

    for t = 1:numel(S.anot_data_all)
        D = S.anot_data_all{t};

        % Align labels to reference order
        allLbl = D.label(:);
        if exist('do_ensure_consistent_labels','file')==2
            D = do_ensure_consistent_labels(D, refLbl);
            allLbl = D.label(:);
        end
        megMask    = startsWith(allLbl,'MEG');
        MEG_labels = allLbl(megMask);

        E = D.trial{1}(megMask, :);          % C×T
        if exist('maggrad_scale','file')==2, E = maggrad_scale(E, MEG_labels); end
        if exist('do_normalize_data','file')==2
            E = do_normalize_data(E,'demean');
        else
            E = E - mean(E,2);
        end

        % Fallback alignment to refLbl if needed
        if ~isequal(MEG_labels(:), refLbl(:))
            m = containers.Map(MEG_labels, 1:numel(MEG_labels));
            idx = zeros(numel(refLbl),1); ok = true;
            for iLab = 1:numel(refLbl)
                if isKey(m, refLbl{iLab}), idx(iLab) = m(refLbl{iLab}); else, ok = false; break; end
            end
            if ~ok
                if VERBOSE, fprintf('[%s] %s: label mismatch, epoch %d skipped.\n', tag, filelist(k).name, t); end
                continue;
            end
            E = E(idx, :);
        end

        cells{end+1,1}  = E;     %#ok<AGROW>
        groups{end+1,1} = gid;   %#ok<AGROW>
        nThis  = nThis + 1;
        totEpochs = totEpochs + 1;
    end

    if VERBOSE
        fprintf('[%s] %s: epochs added = %d (total=%d)\n', tag, filelist(k).name, nThis, totEpochs);
    end
    prog(k,nFiles,totEpochs);
end

    function prog(i,n,tot)
        frac = i/max(n,1);
        if useFT
            ft_progress(frac, '[%s] File %d/%d | total epochs: %d', tag, i, n, tot);
        elseif USEWAIT && ishandle(wb)
            waitbar(frac, wb, sprintf('[%s] Reading %d/%d | epochs: %d', tag, i, n, tot));
        else
            if mod(i, max(1, floor(n/20)))==0  % ~20 updates
                fprintf('[%s] %d/%d files (epochs=%d)\n', tag, i, n, tot);
            end
        end
    end
end


function gid = local_group_id(fullpath)
% token before '_' in filename else parent folder name
[folder, base, ~] = fileparts(fullpath);
[~, parent] = fileparts(folder);
m = regexp(base,'^([^_]+)_','tokens','once');
if ~isempty(m), gid = m{1}; else, gid = parent; end
end


% function trials = build_trials_from_mat(filesS, filesN)
% % Build aligned MEG epoch cells + labels + groups from annotated .mat files.
% % Returns:
% %   trials.data   : {N×1} cell of [C×T] epochs
% %   trials.labels : N×1 categorical ('NoSpike','Spike')
% %   trials.groups : {N×1} group ids (subject/session)
% %   trials.refLbl : C×1 label list (row order used for all epochs)
% 
% % ---- 1) Determine reference label order from the first Spike file ----
% refLbl = [];
% for k = 1:numel(filesS)
%     S = load(filesS(k).name);
%     if isfield(S,'anot_data_all') && ~isempty(S.anot_data_all)
%         D0 = S.anot_data_all{1};
%         allLbl  = D0.label(:);
%         megMask = startsWith(allLbl,'MEG');
%         refLbl  = allLbl(megMask);
%         break
%     end
% end
% assert(~isempty(refLbl),'Could not determine reference MEG labels from spike files.');
% 
% % ---- 2) Read files into epoch cells (Spike first, then NoSpike) ----
% [dataS, groupsS] = read_filelist(filesS, refLbl);
% [dataN, groupsN] = read_filelist(filesN, refLbl);
% 
% % ---- 3) Remove any NaN epochs class-wise (rare but safer) ----
% isBad = @(x) any(isnan(x(:)));
% keepS = ~cellfun(isBad, dataS);
% keepN = ~cellfun(isBad, dataN);
% dataS = dataS(keepS); groupsS = groupsS(keepS);
% dataN = dataN(keepN); groupsN = groupsN(keepN);
% 
% % ---- 4) Concatenate with categorical labels (fixed order) ----
% data   = [dataN; dataS];
% labels = categorical([zeros(numel(dataN),1); ones(numel(dataS),1)], [0 1], {'NoSpike','Spike'});
% groups = [groupsN; groupsS];
% 
% % ---- 5) Sanity checks: same #channels everywhere, matches refLbl ----
% chCounts = cellfun(@(x) size(x,1), data);
% assert(all(chCounts==numel(refLbl)), 'Channel mismatch: not all epochs match refLbl length.');
% 
% % ---- 6) Return struct
% trials = struct('data',{data}, 'labels',labels, 'groups',{groups}, 'refLbl',{refLbl});
% end
% 
% function [cells, groups] = read_filelist(filelist, refLbl)
% % Load all anot_data_all epochs from a list of files, align to refLbl.
% cells  = {};
% groups = {};
% for k = 1:numel(filelist)
%     S = load(filelist(k).name);
%     if ~isfield(S,'anot_data_all') || isempty(S.anot_data_all), continue; end
% 
%     % robust group id (no dependency on extract_group_id)
%     gid = local_group_id(filelist(k).name);
% 
%     for t = 1:numel(S.anot_data_all)
%         D = S.anot_data_all{t};
% 
%         % Align labels to reference order (use helper if available)
%         if exist('do_ensure_consistent_labels','file')==2
%             D = do_ensure_consistent_labels(D, refLbl);
%             allLbl = D.label(:);
%         else
%             allLbl = D.label(:);
%         end
% 
%         % Keep MEG only, scale mags/grads, demean
%         megMask    = startsWith(allLbl,'MEG');
%         MEG_labels = allLbl(megMask);
%         E          = D.trial{1}(megMask, :);          % C×T
%         if exist('maggrad_scale','file')==2, E = maggrad_scale(E, MEG_labels); end
%         if exist('do_normalize_data','file')==2
%             E = do_normalize_data(E,'demean');
%         else
%             E = E - mean(E,2);
%         end
% 
%         % Fallback alignment to refLbl if needed
%         if ~isequal(MEG_labels(:), refLbl(:))
%             m = containers.Map(MEG_labels, 1:numel(MEG_labels));
%             idx = zeros(numel(refLbl),1); ok = true;
%             for iLab = 1:numel(refLbl)
%                 if isKey(m, refLbl{iLab}), idx(iLab) = m(refLbl{iLab}); else, ok = false; break; end
%             end
%             if ~ok, continue; end
%             E = E(idx, :);
%         end
% 
%         cells{end+1,1}  = E;    %#ok<AGROW>
%         groups{end+1,1} = gid;  %#ok<AGROW>
%     end
% end
% end
% 
% function gid = local_group_id(fullpath)
% % token before '_' in filename else parent folder name
% [folder, base, ~] = fileparts(fullpath);
% [~, parent] = fileparts(folder);
% m = regexp(base,'^([^_]+)_','tokens','once');
% if ~isempty(m), gid = m{1}; else, gid = parent; end
% end
