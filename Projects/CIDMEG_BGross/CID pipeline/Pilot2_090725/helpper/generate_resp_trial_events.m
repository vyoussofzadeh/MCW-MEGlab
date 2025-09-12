function generate_resp_trial_events(csvDir, outDir, varargin)
% Pair each trial onset with the FIRST subsequent resp1/resp2 (before next trial,
% within a lag), and write Array-of-times files per run:
%   <run>_trial_resp1.txt   and   <run>_trial_resp2.txt
%
% Name-Value:
%   'TimeUnit'   : 's' (default) or 'ms'
%   'MaxRespLag' : max seconds from trial to response (default 10)
%   'TrialRegex' : default '^trial\d*$'
%   'Resp1Regex' : default '^resp1$'
%   'Resp2Regex' : default '^resp2$'

p = inputParser;
addRequired(p,'csvDir',@(s)ischar(s)||isstring(s));
addRequired(p,'outDir',@(s)ischar(s)||isstring(s));
addParameter(p,'TimeUnit','s',@(s)any(strcmpi(s,{'s','ms'})));
addParameter(p,'MaxRespLag',10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TrialRegex','^trial\d*$',@ischar);
addParameter(p,'Resp1Regex','^resp1$',@ischar);
addParameter(p,'Resp2Regex','^resp2$',@ischar);
parse(p,csvDir,outDir,varargin{:});
opt = p.Results;

csvDir = char(csvDir); outDir = char(outDir);
if ~exist(outDir,'dir'), mkdir(outDir); end

D = dir(fullfile(csvDir, '*.csv'));
if isempty(D), error('No CSV files found in %s', csvDir); end

for i = 1:numel(D)
    inFile = fullfile(D(i).folder, D(i).name);
    [~, base] = fileparts(inFile);

    % --- Read first two columns robustly (label,time) ---
    [lab, t] = read_label_time_compat(inFile);
    if strcmpi(opt.TimeUnit,'ms'), t = t/1000; end

    % Clean + sort
    good = lab~="" & ~isnan(t);
    lab = lower(strtrim(lab(good))); t = t(good);
    [t, idx] = sort(t); lab = lab(idx);

    % Classify
    isTrial = ~cellfun('isempty', regexp(cellstr(lab), opt.TrialRegex, 'once'));
    isR1    = ~cellfun('isempty', regexp(cellstr(lab), opt.Resp1Regex, 'once'));
    isR2    = ~cellfun('isempty', regexp(cellstr(lab), opt.Resp2Regex, 'once'));

    t_trial = t(isTrial);
    t_r1    = t(isR1);
    t_r2    = t(isR2);

    % Pair trials with first subsequent response (before next trial, within lag)
    t_resp1_trials = [];
    t_resp2_trials = [];

    if ~isempty(t_trial)
        nextTrial = [t_trial(2:end); +inf];
        for k = 1:numel(t_trial)
            t0 = t_trial(k); tNext = nextTrial(k);

            c1 = t_r1(t_r1 > t0 & t_r1 < tNext & (t_r1 - t0) <= opt.MaxRespLag);
            c2 = t_r2(t_r2 > t0 & t_r2 < tNext & (t_r2 - t0) <= opt.MaxRespLag);

            % Lazy selection (no c1(1)/c2(1) unless non-empty)
            if ~isempty(c1)
                t1 = c1(1);
            else
                t1 = inf;
            end
            if ~isempty(c2)
                t2 = c2(1);
            else
                t2 = inf;
            end

            if t1 < t2
                t_resp1_trials(end+1,1) = t0; %#ok<AGROW>
            elseif t2 < t1
                t_resp2_trials(end+1,1) = t0; %#ok<AGROW>
            else
                % no response within window => timeout (ignored)
            end
        end
    end

    % --- Write Array-of-times files ---
    f1 = fullfile(outDir, sprintf('%s_trial_resp1.txt', base));
    f2 = fullfile(outDir, sprintf('%s_trial_resp2.txt', base));
    write_times(f1, t_resp1_trials);
    write_times(f2, t_resp2_trials);

    fprintf('[%s]\n  trials=%d  resp1-trials=%d  resp2-trials=%d  (timeouts=%d)\n', ...
        base, numel(t_trial), numel(t_resp1_trials), numel(t_resp2_trials), ...
        numel(t_trial) - numel(t_resp1_trials) - numel(t_resp2_trials));
end
end

% ---------- helpers ----------
function [lab, t] = read_label_time_compat(inFile)
    % Try modern readtable (2023b supports VariableNamingRule, but simple call is fine)
    try
        T = readtable(inFile);
        if width(T) >= 2
            lab = string(T{:,1});
            t   = double(T{:,2});
            return;
        end
    catch
    end
    % Fallback: textscan (comma/semicolon/tab/space)
    fid = fopen(inFile,'r'); assert(fid>0, 'Cannot open %s', inFile);
    C = textscan(fid, '%s %f', 'Delimiter', {',',';','\t',' '}, ...
                 'MultipleDelimsAsOne', true, 'HeaderLines', 0);
    fclose(fid);
    lab = string(C{1});
    t   = double(C{2});
end

function write_times(fname, t)
    fid = fopen(char(fname),'w'); assert(fid>0,'Cannot write %s', fname);
    for ii = 1:numel(t), fprintf(fid,'%.6f\n', t(ii)); end
    fclose(fid);
end