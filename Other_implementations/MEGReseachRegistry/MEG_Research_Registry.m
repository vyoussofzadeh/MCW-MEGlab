%% MEG scan list + dates (flat list, unique dates, counts, exports)
% Author: VY
% Last updated: 2025-10-28
% NOTE: Older-MATLAB compatible (no groupcounts)

%% --- Paths
clear; clc;
datadir = '/MEG_data/Research_studies/MEG_Reseach_records';
cd(datadir);

din = fullfile(datadir, '/Imported/2025/Research_scans_anncal_10282025.xlsx');

%% --- Select year
Y_var = input('enter year (e.g., 2025): ', 's');   % <---- change this when needed
Y_var = str2double(Y_var);

%% --- Read and normalize table
opts = detectImportOptions(din);

% Ensure Start/End come in as text so we can clean robustly
mustHave = {'Start','End'};
for v = mustHave
    if ismember(v{1}, opts.VariableNames)
        opts = setvartype(opts, v{1}, 'char');
    end
end

data = readtable(din, opts);

% Find Start/End columns case-insensitively
vn = string(data.Properties.VariableNames);
startVar = vn(ismember(lower(vn), "start"));
endVar   = vn(ismember(lower(vn), "end"));
if isempty(startVar) || isempty(endVar)
    error('Could not find Start/End columns (case-insensitive).');
end
startVar = char(startVar(1));
endVar   = char(endVar(1));

% --- Normalize strings aggressively (helps 2024 rows)
data.(startVar) = strtrim(data.(startVar));
data.(endVar)   = strtrim(data.(endVar));

% Remove leading weekday names + optional comma (e.g., "Tue, 10/03/24 ...")
data.(startVar) = regexprep(data.(startVar), '^\s*\w{3,9},?\s+', '');
data.(endVar)   = regexprep(data.(endVar),   '^\s*\w{3,9},?\s+', '');

% Remove tokens like " at ", timezone labels (CT|CST|CDT), UTC/GMT
killTokens = {'\sat\s','\sCT\b','\sCST\b','\sCDT\b','\sUTC[+-]?\d*','\sGMT\b'};
for k = 1:numel(killTokens)
    data.(startVar) = regexprep(data.(startVar), killTokens{k}, ' ', 'ignorecase');
    data.(endVar)   = regexprep(data.(endVar),   killTokens{k}, ' ', 'ignorecase');
end

% Collapse multiple spaces
data.(startVar) = regexprep(data.(startVar), '\s+', ' ');
data.(endVar)   = regexprep(data.(endVar),   '\s+', ' ');

% Common datetime formats we’ll try (covers 12h/24h, seconds, month names)
fmts = { ...
    'MM/dd/yyyy h:mm a',    'M/d/yyyy h:mm a', ...
    'MM/dd/yy h:mm a',      'M/d/yy h:mm a', ...
    'MM/dd/yyyy HH:mm',     'M/d/yyyy HH:mm', ...
    'MM/dd/yy HH:mm',       'M/d/yy HH:mm', ...
    'MM/dd/yyyy h:mm:ss a', 'M/d/yyyy h:mm:ss a', ...
    'MM/dd/yy h:mm:ss a',   'M/d/yy h:mm:ss a', ...
    'MM/dd/yyyy',           'M/d/yyyy', ...
    'MMM d, yyyy h:mm a',   'MMMM d, yyyy h:mm a', ...
    'MMM d, yyyy',          'MMMM d, yyyy' ...
};

% Parse to datetime with fallbacks (handles Excel serials too)
data.Start = parseDatetimeWithFallback(data.(startVar), fmts, 'en_US');
data.End   = parseDatetimeWithFallback(data.(endVar),   fmts, 'en_US');

% Drop rows with unparsed times
bad = isnat(data.Start) | isnat(data.End);
if any(bad)
    warning('Dropping %d rows with unparsed Start/End times.', sum(bad));
    data(bad,:) = [];
end

% Duration
data.Duration = data.End - data.Start;

% Subject must exist
if ~ismember('Subject', data.Properties.VariableNames)
    error('Expected a "Subject" column in the file.');
end

%% --- Date-only & slice selected year
data.ScanDate = dateshift(data.Start, 'start', 'day');

yrMask = data.ScanDate >= datetime(Y_var,1,1) & data.ScanDate < datetime(Y_var+1,1,1);
dataY  = data(yrMask,:);

% Add helpful fields
dataY.Day     = categorical(cellstr(datestr(dataY.ScanDate, 'yyyy-mm-dd (ddd)')));
dataY.StartTT = timeofday(dataY.Start);
dataY.EndTT   = timeofday(dataY.End);
dataY.Hours   = hours(dataY.Duration);

%% --- Flat list
flatList = sortrows(dataY(:, {'Subject','ScanDate','Day','Start','End','StartTT','EndTT','Hours'}), ...
                    {'Subject','ScanDate','Start'});

disp(['--- ' num2str(Y_var) ' scans (flat list) ---']);
disp(flatList);

%% --- Unique list: one row per Subject × ScanDate
[uniqPairs, ia] = unique(flatList(:, {'Subject','ScanDate','Day'}), 'rows', 'stable');
uniqueDates = sortrows(uniqPairs, {'Subject','ScanDate'});

disp(['--- ' num2str(Y_var) ' unique scan dates per project ---']);
disp(uniqueDates);

%% --- Per-project daily counts (older-MATLAB style)
tmpCounts = varfun(@numel, flatList, ...
    'InputVariables','ScanDate', ...
    'GroupingVariables',{'Subject','ScanDate'});

% Find the count column produced by varfun
vn = tmpCounts.Properties.VariableNames;
idx = find(strcmp(vn,'numel_ScanDate'), 1);
if isempty(idx), idx = find(strncmp(vn,'numel_',6), 1, 'last'); end
if isempty(idx), error('Could not locate the count column from varfun.'); end

tmpCounts.Properties.VariableNames{idx} = 'NumBookings';
perDayCounts = sortrows(tmpCounts, {'Subject','ScanDate'});

disp(['--- ' num2str(Y_var) ' counts per project/day ---']);
disp(perDayCounts);

%% --- Per-project totals across the year (older-MATLAB style)
tmpTotals = varfun(@numel, flatList, ...
    'InputVariables','ScanDate', ...
    'GroupingVariables','Subject');

vn  = tmpTotals.Properties.VariableNames;
idx = find(strcmp(vn,'numel_ScanDate'), 1);
if isempty(idx), idx = find(strncmp(vn,'numel_',6), 1, 'last'); end
if isempty(idx), error('Could not locate the count column for totals.'); end

tmpTotals.Properties.VariableNames{idx} = ['NumBookings' num2str(Y_var)];
try
    perProjectTotals = sortrows(tmpTotals, tmpTotals.Properties.VariableNames{idx}, 'descend');
catch
    perProjectTotals = sortrows(tmpTotals, idx, 'descend');
end

disp(['--- ' num2str(Y_var) ' totals per project ---']);
disp(perProjectTotals);

%% --- Quarter tallies (OCT–DEC of previous year, then Q1/Q2 of selected year)
quarters = {
    ['OCT–DEC ' num2str(Y_var-1)], datetime(Y_var-1,10,1), datetime(Y_var,1,1);
    ['JAN–MAR ' num2str(Y_var)],   datetime(Y_var,1,1),    datetime(Y_var,4,1);
    ['APR–JUN ' num2str(Y_var)],   datetime(Y_var,4,1),    datetime(Y_var,7,1);
};

fprintf('\n=== Quarter summaries (relative to %d) ===\n', Y_var);
for q = 1:size(quarters,1)
    label = quarters{q,1};
    t0    = quarters{q,2};
    t1    = quarters{q,3};

    m = data.Start >= t0 & data.Start < t1;         % across raw data to include OCT–DEC (Y-1)
    nScans = sum(m);
    totHrs = sum(hours(data.Duration(m)));
    nDays  = numel(unique(dateshift(data.Start(m), 'start', 'day')));

    fprintf('%s = %d scans, %d unique days, %.1f hours\n', label, nScans, nDays, totHrs);
end

%% --- Exports (XLSX with multiple sheets + CSVs), named by selected year
outDir  = fullfile(datadir, 'Shared');
if ~exist(outDir, 'dir'), mkdir(outDir); end

outXLSX = fullfile(outDir, sprintf('MEG_scans_%d_list.xlsx', Y_var));
outCSV1 = fullfile(outDir, sprintf('MEG_scans_%d_list_flat.csv', Y_var));
outCSV2 = fullfile(outDir, sprintf('MEG_scans_%d_unique_dates.csv', Y_var));
outCSV3 = fullfile(outDir, sprintf('MEG_scans_%d_counts_per_day.csv', Y_var));
outCSV4 = fullfile(outDir, sprintf('MEG_scans_%d_totals.csv', Y_var));

try
    if exist(outXLSX, 'file'); delete(outXLSX); end
catch
    % ignore delete failure
end

writetable(flatList,         outXLSX, 'Sheet', 'flat');
writetable(uniqueDates,      outXLSX, 'Sheet', 'unique_dates');
writetable(perDayCounts,     outXLSX, 'Sheet', 'counts_per_day');
writetable(perProjectTotals, outXLSX, 'Sheet', 'totals');

writetable(flatList,         outCSV1);
writetable(uniqueDates,      outCSV2);
writetable(perDayCounts,     outCSV3);
writetable(perProjectTotals, outCSV4);

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n  %s\n', outXLSX, outCSV1, outCSV2, outCSV3, outCSV4);

%% --- Helper: robust datetime parser with fallbacks (serials, formats, free-parse)
function dt = parseDatetimeWithFallback(txtCol, fmtList, locale)
    % Accept cellstr/string/char/numeric and convert robustly.
    if isstring(txtCol), txtCol = cellstr(txtCol); end
    if isnumeric(txtCol)
        dt = datetime(txtCol, 'ConvertFrom','excel');  % Excel serials
        return
    end
    if iscell(txtCol)
        c = txtCol;
        for i = 1:numel(c)
            if isnumeric(c{i}) && ~isempty(c{i})
                c{i} = num2str(c{i}, 16);
            elseif isstring(c{i})
                c{i} = char(c{i});
            elseif ischar(c{i})
                % ok
            elseif isempty(c{i})
                c{i} = '';
            else
                c{i} = char(string(c{i}));
            end
        end
        txtCol = c;
    elseif ischar(txtCol)
        txtCol = cellstr(txtCol);
    end

    dt = NaT(size(txtCol));
    mask = isnat(dt);

    % Try Excel serials hiding as strings (e.g., '45231.5')
    isSerial = false(size(txtCol));
    for i = 1:numel(txtCol)
        t = strtrim(txtCol{i});
        if ~isempty(t) && ~isnan(str2double(t)) && all(regexp(t, '^[0-9]+(\.[0-9]+)?$'))
            isSerial(i) = true;
        end
    end
    if any(isSerial)
        serialVals = cellfun(@str2double, txtCol(isSerial));
        dt(isSerial) = datetime(serialVals, 'ConvertFrom','excel');
        mask = isnat(dt);
    end

    % Try the provided formats in order
    for k = 1:numel(fmtList)
        if ~any(mask), break; end
        try
            trial = datetime(txtCol(mask), 'InputFormat', fmtList{k}, 'Locale', locale);
        catch
            trial = NaT(sum(mask),1);
        end
        dt(mask) = trial;
        mask = isnat(dt);
    end

    % Last-chance free parse for month-name strings
    if any(mask)
        try
            trial = datetime(txtCol(mask), 'Locale', locale);
            dt(mask) = trial;
        catch
            % leave as NaT
        end
    end
end
