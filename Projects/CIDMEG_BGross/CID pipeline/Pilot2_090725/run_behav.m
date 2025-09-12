% tag_trials_from_behavior.m
% Tag Brainstorm trials using a behavioral table (run, trial, key, rt_category)
% Run this in the directory with data_*_trial*.mat files.

clear; clc;

% ---- CONFIG ----
behavFile = 'CIDMEG004_behavioral_data.csv';     % <- put your table filename here
markTimeoutAsBad = true;              % add " [bad]" to comments for timeouts
fastCode = 1; slowCode = 2; timeoutCode = 3;

% ---- LOAD TABLE ----
T = readtable(behavFile, 'FileType','text', 'Delimiter', detectImportOptions(behavFile).Delimiter);

% Normalize columns (handle 'timeout' string in key)
if ~isnumeric(T.key)
    % Map 'timeout' -> NaN (or 0), numeric strings -> numbers
    keyNum = nan(height(T),1);
    for i = 1:height(T)
        if ismember(lower(string(T.key(i))), ["timeout","time_out","na",""])
            keyNum(i) = NaN;
        else
            keyNum(i) = str2double(string(T.key(i)));
        end
    end
    T.key = keyNum;
end

% Sanity check
reqCols = ["run","trial","rt_category","key"];
missingCols = setdiff(reqCols, string(T.Properties.VariableNames));
if ~isempty(missingCols)
    error('Missing required columns in table: %s', strjoin(missingCols, ', '));
end

% ---- HELPERS ----
mkFile = @(r,t) sprintf('data_%d_trial%03d.mat', r, t);
rtTag  = @(c)  ternary(c==fastCode, 'RT=fast', ternary(c==slowCode,'RT=slow','RT=timeout'));
keyTag = @(k)  ( ~isnan(k) .* k + isnan(k).*0 ); % numeric or 0 for NaN
hasStr = @(s,sub) ~isempty(strfind(s, sub)); %#ok<STREMP>

% ---- TAG TRIALS ----
nFound=0; nMissing=0; nFast=0; nSlow=0; nTimeout=0;

for i = 1:height(T)
    runNo   = T.run(i);
    trialNo = T.trial(i);
    rtCat   = T.rt_category(i);
    keyVal  = T.key(i);
    stim    = string(T.stimulus(i));

    fn = mkFile(runNo, trialNo);
    if ~isfile(fn)
        fprintf('MISS: %s\n', fn);
        nMissing = nMissing + 1;
        continue;
    end

    S = load(fn); % DataMat fields: F, ChannelFlag, Comment, Time, ...
    if ~isfield(S,'Comment'), S.Comment = ''; end

    % Compose tags
    tags = sprintf(' |Run=%d| |Trial=%03d| |%s|', runNo, trialNo, rtTag(rtCat));
    if ~isnan(keyVal)
        tags = sprintf('%s |Key=%d|', tags, keyVal);
    else
        tags = sprintf('%s |Key=NA|', tags);
    end
    % Optional stimulus base tag (without _1/_2/_3)
    baseStim = regexprep(stim, '_\d+\.png$', '');
    if ~isempty(baseStim)
        tags = sprintf('%s |Stim=%s|', tags, baseStim);
    end

    % Append tags only if not already present
    newComment = S.Comment;
    for token = regexp(strtrim(tags), '\|[^|]+\|', 'match')
        tok = token{1};
        if ~contains(newComment, tok), newComment = [newComment ' ' tok]; end %#ok<AGROW>
    end

    % Mark timeout as bad (optional)
    if markTimeoutAsBad && (rtCat == timeoutCode)
        if ~contains(newComment, '[bad]'), newComment = [newComment ' [bad]']; end
    end

    % Update counts
    switch rtCat
        case fastCode,    nFast = nFast + 1;
        case slowCode,    nSlow = nSlow + 1;
        otherwise,        nTimeout = nTimeout + 1;
    end

    % Append to History if present
    if isfield(S,'History') && iscell(S.History)
        try
            S.History(end+1,1:3) = {datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
                                    'tag_trials_from_behavior', strtrim(tags)};
        catch
            % ignore
        end
    end

    % Save
    S.Comment = strtrim(newComment);
    save(fn, '-struct', 'S');

    fprintf('TAG : %s  -> %s\n', fn, S.Comment);
    nFound = nFound + 1;
end

fprintf('\nSummary: Tagged=%d, Missing=%d | Fast=%d, Slow=%d, Timeout=%d\n', ...
        nFound, nMissing, nFast, nSlow, nTimeout);


% Tiny ternary
function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end



