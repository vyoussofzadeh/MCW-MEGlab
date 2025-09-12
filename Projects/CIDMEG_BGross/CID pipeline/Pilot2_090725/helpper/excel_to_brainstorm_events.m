function excel_to_brainstorm_events(inputPath, varargin)
% excel_to_brainstorm_events  Convert Excel/CSV event tables to Brainstorm text.
%
% INPUT
%   inputPath : file (.xlsx/.xls/.csv/.txt) or folder.
%
% Name-Value (optional)
%   'Format'    : 'array' | 'trg' | 'both'    (default 'both')
%   'TimeUnit'  : 's' | 'ms'                  (default 's')
%   'Fs'        : sampling rate in Hz for .trg sample column (default [], write 0)
%   'LabelCol'  : header name of label/type column (default: auto)
%   'TimeCol'   : header name of time/onset column  (default: auto)
%   'LabelIdx'  : 1-based column index for label    (default: auto)
%   'TimeIdx'   : 1-based column index for time     (default: auto)
%   'OutputDir' : output root directory (default: <inDir>/bst_events/<filebase>)
%
% Brainstorm import:
%   - Array-of-times: Record tab ? File ? Add events from file ? "Array of times"
%   - .trg:           Record tab ? File ? Add events from file ? "ANT EEProbe (.trg)"

p = inputParser;
addRequired(p, 'inputPath', @(s)ischar(s)||isstring(s));
addParameter(p, 'Format',   'both', @(s)any(strcmpi(s,{'array','trg','both'})));
addParameter(p, 'TimeUnit', 's',    @(s)any(strcmpi(s,{'s','ms'})));
addParameter(p, 'Fs',       [],     @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
addParameter(p, 'LabelCol', '',     @(s)ischar(s)||isstring(s));
addParameter(p, 'TimeCol',  '',     @(s)ischar(s)||isstring(s));
addParameter(p, 'LabelIdx', [],     @(x) isempty(x) || (isscalar(x)&&x>=1));
addParameter(p, 'TimeIdx',  [],     @(x) isempty(x) || (isscalar(x)&&x>=1));
addParameter(p, 'OutputDir','',     @(s)ischar(s)||isstring(s));
parse(p, inputPath, varargin{:});
opt = p.Results;

fmtArray  = any(strcmpi(opt.Format, {'array','both'}));
fmtTrg    = any(strcmpi(opt.Format, {'trg','both'}));
unitScale = strcmpi(opt.TimeUnit,'ms')*(1/1000) + strcmpi(opt.TimeUnit,'s')*1;

% --- Collect files
if isfolder(opt.inputPath)
    D = dir(fullfile(opt.inputPath, '*'));
    fileList = {};
    for k = 1:numel(D)
        if D(k).isdir, continue; end
        [~,~,e] = fileparts(D(k).name);
        if any(strcmpi(e,{'.xlsx','.xls','.csv','.txt'}))
            fileList{end+1} = fullfile(D(k).folder, D(k).name); %#ok<AGROW>
        end
    end
else
    fileList = {char(opt.inputPath)};
end
if isempty(fileList), error('No input files found.'); end
fprintf('Found %d file(s).\n', numel(fileList));

for i = 1:numel(fileList)
    inFile = fileList{i};
    [inDir, base] = fileparts(inFile);

    % Output dir
    if ~isempty(opt.OutputDir)
        outRoot = char(opt.OutputDir);
    else
        outRoot = fullfile(inDir, 'bst_events');
    end
    outDir = fullfile(outRoot, base);
    if ~exist(outDir,'dir'), mkdir(outDir); end

    % Read table (let MATLAB detect delimiter/headers; robust to odd headers)
    try
        T = readtable(inFile, 'VariableNamingRule','preserve');
    catch
        T = readtable(inFile, 'ReadVariableNames', false);  % last resort
    end
    varNames = string(T.Properties.VariableNames);
    nvars    = numel(varNames);

    % ----- Decide label/time columns -----
    [labelVar, timeVar] = pickColumns(T, varNames, opt);

    % Pull data
    if ~isempty(opt.LabelIdx)
        labelsRaw = T{:, opt.LabelIdx};
    else
        labelsRaw = T.(labelVar);
    end
    if ~isempty(opt.TimeIdx)
        timesRaw = T{:, opt.TimeIdx};
    else
        timesRaw = T.(timeVar);
    end

    labels = toStrCol(labelsRaw);
    times  = toDoubleCol(timesRaw);

    % Apply time unit
    times = times .* unitScale;

    % Drop invalids
    good = ~isnan(times) & labels ~= "";
    if ~all(good)
        fprintf('  %s: dropping %d invalid row(s).\n', base, sum(~good));
    end
    labels = labels(good);
    times  = times(good);

    % Sort
    [times, idx] = sort(times);
    labels = labels(idx);

    if isempty(times)
        fprintf('  %s: no valid events, skipping.\n', base);
        continue;
    end

    % ----- Write Array-of-times (per label)
    if fmtArray
        labU  = unique(labels,'stable');
        labDir = fullfile(outDir,'array_of_times');
        if ~exist(labDir,'dir'), mkdir(labDir); end
        for j = 1:numel(labU)
            lab = labU(j);
            tLab = times(labels==lab);

            % fn = fullfile(labDir, [sanitizeForFile(lab) '.txt']);
            % fid = fopen(fn,'w'); assert(fid>0, 'Cannot write %s', fn);
            % fprintf(fid, '%.6f\n', tLab);
            % fclose(fid);

            labName = sanitizeForFile(lab);                  % char
            fn = fullfile(labDir, [labName '.txt']);        % char
            fid = fopen(fn,'w'); assert(fid>0, 'Cannot write %s', fn);
            fprintf(fid, '%.6f\n', tLab);
            fclose(fid);

        end
        fprintf('  %s: wrote Array-of-times files (%d labels).\n', base, numel(labU));
    end

    % ----- Write combined .trg
    if fmtTrg
        trgDir  = fullfile(outDir,'trg');
        if ~exist(trgDir,'dir'), mkdir(trgDir); end
        trgFile = fullfile(trgDir,'events_combined.trg');
        fid = fopen(trgFile,'w'); assert(fid>0, 'Cannot write %s', trgFile);
        fprintf(fid, 'latency sample name\n');
        if isempty(opt.Fs)
            for k = 1:numel(times)
                fprintf(fid, '%.6f %d %s\n', times(k), 0, labels(k));
            end
        else
            samp = round(times .* opt.Fs);
            for k = 1:numel(times)
                fprintf(fid, '%.6f %d %s\n', times(k), samp(k), labels(k));
            end
        end
        fclose(fid);
        fprintf('  %s: wrote .trg file.\n', base);
    end
end
fprintf('Done.\n');

% --------- helpers ---------
function [labelVar, timeVar] = pickColumns(T, varNames, opt)
    % 1) If user provided indices, we can skip header lookup
    if ~isempty(opt.LabelIdx) && ~isempty(opt.TimeIdx)
        labelVar = varNames(max(1,min(end,opt.LabelIdx)));
        timeVar  = varNames(max(1,min(end,opt.TimeIdx)));
        return;
    end

    % 2) If user provided header names
    labelVar = string(opt.LabelCol);
    timeVar  = string(opt.TimeCol);

    % 3) Auto-detect by common names
    if labelVar=="" || ~any(varNames==labelVar)
        labelVar = firstMatch(varNames, ["label","Label","event","Event","type","Type","name","Name"]);
    end
    if timeVar=="" || ~any(varNames==timeVar)
        timeVar = firstMatch(varNames, ["time","Time","onset","Onset","latency","Latency","t","ms","Msec","msec"]);
    end

    % 4) Final fallback: first two columns by position
    if (labelVar=="" || ~any(varNames==labelVar)) || (timeVar=="" || ~any(varNames==timeVar))
        if numel(varNames) < 2
            error('Need at least two columns (label,time).');
        end
        labelVar = varNames(1);
        timeVar  = varNames(2);
    end
end

function out = firstMatch(list, keys)
    out = "";
    for kk = 1:numel(keys)
        if any(list==keys(kk)), out = keys(kk); return; end
    end
end

function s = sanitizeForFile(s)
    % s = regexprep(string(s), '[^a-zA-Z0-9_\-]', '_');
    s = char(regexprep(string(s), '[^a-zA-Z0-9_\-]', '_'));
end

function s = toStrCol(x)
    if isstring(x), s = x; return; end
    if iscellstr(x), s = string(x); return; end
    if iscell(x)
        s = strings(size(x));
        for ii=1:numel(x), s(ii) = string(x{ii}); end
        return;
    end
    if ischar(x), s = string(x); return; end
    s = string(x);  % numeric ? string
end

function v = toDoubleCol(x)
    if isnumeric(x), v = double(x); return; end
    if iscell(x)
        v = nan(size(x));
        for ii=1:numel(x), v(ii) = str2double(string(x{ii})); end
        return;
    end
    v = str2double(string(x));
end
end
