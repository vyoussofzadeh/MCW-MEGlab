function varargout = process_average_timeinterval( varargin )
% PROCESS_AVERAGE_TIMEINTERVAL: For each file in input, compute the mean (or
% the variance, etc.) over multiple time intervals within the time window.

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Average time interval';
sProcess.FileTag     = @GetFileTag;
sProcess.Category    = 'Filter';
sProcess.SubGroup    = 'Average';
sProcess.Index       = 304; % Ensure unique index
sProcess.Description = '';
% Input/output types
sProcess.InputTypes  = {'data', 'results', 'timefreq', 'matrix'};
sProcess.OutputTypes = {'data', 'results', 'timefreq', 'matrix'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

% New options for interval configuration
sProcess.options.strt.Comment = 'Start time (s):';
sProcess.options.strt.Type    = 'value';
sProcess.options.strt.Value   = {0, 's', 2};

sProcess.options.spt.Comment = 'Span time (s):';
sProcess.options.spt.Type    = 'value';
sProcess.options.spt.Value   = {2, 's', 2};

sProcess.options.overlap.Comment = 'Overlap (s):';
sProcess.options.overlap.Type    = 'value';
sProcess.options.overlap.Value   = {0.01, 's', 3};

sProcess.options.linterval.Comment = 'Interval length (s):';
sProcess.options.linterval.Type    = 'value';
sProcess.options.linterval.Value   = {0.5, 's', 2};

% Inherit other options from original process_average_time
% ...
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
% Get time window
Comment = [upper(GetFileTag(sProcess)), ': [', process_extract_time('GetTimeString', sProcess), ']'];
% Absolute values
if isfield(sProcess.options, 'source_abs') && sProcess.options.source_abs.Value
    Comment = [Comment, ', abs'];
end
end

%% ===== GET FILE TAG =====
function fileTag = GetFileTag(sProcess)
% Old version of the process: option isstd={0,1}
if isfield(sProcess.options, 'isstd') && ~isempty(sProcess.options.isstd) && sProcess.options.isstd.Value
    fileTag = 'std';
    % New version of the process: avg_fun={'mean','std','rms','median'}
elseif isfield(sProcess.options, 'avg_func') && ~isempty(sProcess.options.avg_func) && ~isempty(sProcess.options.avg_func.Value)
    fileTag = sProcess.options.avg_func.Value;
else
    fileTag = 'mean';
end
end

%% ===== RUN =====
function sInput = Run(sProcess, sInput) %#ok<DEFNU>
% Calculate time intervals based on process options
cfg_main.strt = sProcess.options.strt.Value{1};
cfg_main.spt = sProcess.options.spt.Value{1};
cfg_main.overlap = sProcess.options.overlap.Value{1};
cfg_main.linterval = sProcess.options.linterval.Value{1};

% Calculate the windows
wi = do_time_intervals(cfg_main);

% Loop through each window and calculate average (or other metric)
for wIdx = 1:size(wi, 1)
    iTime = panel_time('GetTimeIndices', sInput.TimeVector, wi(wIdx,:));
    if isempty(iTime)
        bst_report('Error', sProcess, [], 'Invalid time definition.');
        sInput = [];
        return;
    end
    
    % Apply function for each interval
    switch GetFileTag(sProcess)
        case 'mean'
            tempData = mean(sInput.A(:,iTime,:), 2);
        case 'rms'
            tempData = sqrt(sum(sInput.A(:,iTime,:).^2, 2) / length(iTime));
        case 'std'
            tempData = sqrt(var(sInput.A(:,iTime,:), 0, 2));
        case 'median'
            tempData = median(sInput.A(:,iTime,:), 2);
    end
    
    % Assuming tempData is [numVertices x 1 x numConditions], remove unnecessary dimensions
    tempData = squeeze(tempData);
    
    % Assuming DataMat.Time is sorted and wi intervals do not overlap
    % Compute start and end indices for all intervals in one pass
    [~, start_indices] = arrayfun(@(x) min(abs(sInput.TimeVector - x)), wi(wIdx,1));
    [~, end_indices] = arrayfun(@(x) min(abs(sInput.TimeVector - x)), wi(wIdx,2));
    
    % Populate the matrix with computed averages
    if start_indices < end_indices
        
        avgSourceData(:, start_indices:end_indices) = repmat(tempData, 1, end_indices-start_indices+1);
    end
end

% Update sInput.A with interpolated data
sInput.A = avgSourceData;

% Update TimeVector if necessary
[~, start_indices] = arrayfun(@(x) min(abs(sInput.TimeVector - x)), wi(1,1));
[~, end_indices] = arrayfun(@(x) min(abs(sInput.TimeVector - x)), wi(end,2));
sInput.TimeVector = sInput.TimeVector(start_indices:end_indices); % Assume computeMidpoints is a function to calculate midpoints
sInput.A = sInput.A(:,start_indices:end_indices);

% Update CommentTag to reflect interval processing with interpolation
sInput.CommentTag = [GetFileTag(sProcess) ' with interpolation (' process_extract_time('GetTimeString',sProcess,sInput) ')'];

% Clear unnecessary fields
if isfield(sInput, 'Std') && ~isempty(sInput.Std)
    sInput.Std = [];
end
if isfield(sInput, 'TFmask') && ~isempty(sInput.TFmask)
    sInput.TFmask = [];
end
end


%% ===== DO TIME INTERVALS =====
% Incorporate the do_time_intervals function here
function [wi]  = do_time_intervals(cfg_main)
strt = cfg_main.strt; % strt = 0 sec.
spt = cfg_main.spt; % spt = 2 sec.
overlap = cfg_main.overlap; % overlap = 0.01;
linterval = cfg_main.linterval; % interval length

wi = []; w1 = strt; l = linterval; ov = overlap; j=1;
while w1+l <= spt
    wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
end
end
