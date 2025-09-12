function T = build_trial_response_table(csvDir, outCsv, varargin)
% build_trial_response_table
%   Create a single table mapping (run, trial_idx) -> resp class ('resp1','resp2','timeout').
%   INPUTS:
%     csvDir : folder with per-run CSVs (2 columns: label,time)
%     outCsv : path to save the combined table as CSV (optional; '' to skip)
%   NAME-VALUE:
%     'TimeUnit'   : 's' (default) or 'ms'
%     'MaxRespLag' : max seconds after trial to accept a response (default 10)
%     'TrialRegex' : default '^trial\d*$'
%     'Resp1Regex' : default '^resp1$'
%     'Resp2Regex' : default '^resp2$'
%
%   OUTPUT:
%     T : table(run, trial_idx, resp), plus metadata printed to console.

p = inputParser;
addRequired(p,'csvDir',@(s)ischar(s)||isstring(s));
addRequired(p,'outCsv',@(s)ischar(s)||isstring(s)||isempty(s));
addParameter(p,'TimeUnit','s',@(s)any(strcmpi(s,{'s','ms'})));
addParameter(p,'MaxRespLag',10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TrialRegex','^trial\d*$',@ischar);
addParameter(p,'Resp1Regex','^resp1$',@ischar);
addParameter(p,'Resp2Regex','^resp2$',@ischar);
parse(p,csvDir,outCsv,varargin{:});
opt = p.Results;

D = dir(fullfile(char(opt.csvDir), '*.csv'));
if isempty(D), error('No CSV files found in %s', opt.csvDir); end

rows_run = {}; rows_trial = []; rows_resp = {};

fprintf('Building (run, trial_idx) -> resp table from %d CSV(s)...\n', numel(D));
for i = 1:numel(D)
    inFile = fullfile(D(i).folder, D(i).name);
    [~, runBase] = fileparts(inFile);

    % Read first two columns robustly
    [lab, t] = i_read_label_time(inFile);
    if strcmpi(opt.TimeUnit,'ms'), t = t/1000; end

    good = lab~="" & ~isnan(t);
    lab = lower(strtrim(lab(good))); t = t(good);
    [t, idx] = sort(t); lab = lab(idx);

    isTrial = ~cellfun('isempty', regexp(cellstr(lab), opt.TrialRegex, 'once'));
    isR1    = ~cellfun('isempty', regexp(cellstr(lab), opt.Resp1Regex, 'once'));
    isR2    = ~cellfun('isempty', regexp(cellstr(lab), opt.Resp2Regex, 'once'));

    t_trial = t(isTrial);
    t_r1    = t(isR1);
    t_r2    = t(isR2);

    nextTrial = [t_trial(2:end); +inf];
    trial_idx = (1:numel(t_trial)).';

    resp = strings(size(trial_idx));
    for k = 1:numel(t_trial)
        t0 = t_trial(k); tNext = nextTrial(k);
        c1 = t_r1(t_r1 > t0 & t_r1 < tNext & (t_r1 - t0) <= opt.MaxRespLag);
        c2 = t_r2(t_r2 > t0 & t_r2 < tNext & (t_r2 - t0) <= opt.MaxRespLag);
        t1 = ifelse_first(c1, inf);
        t2 = ifelse_first(c2, inf);
        if t1 < t2
            resp(k) = "resp1";
        elseif t2 < t1
            resp(k) = "resp2";
        else
            resp(k) = "timeout";
        end
    end

    rows_run   = [rows_run; repmat({runBase}, numel(trial_idx), 1)]; %#ok<AGROW>
    rows_trial = [rows_trial; trial_idx]; %#ok<AGROW>
    rows_resp  = [rows_resp;  resp]; %#ok<AGROW>

    fprintf('  %s: trials=%d  resp1=%d  resp2=%d  timeout=%d\n', ...
        runBase, numel(trial_idx), sum(resp=="resp1"), sum(resp=="resp2"), sum(resp=="timeout"));
end

T = table(string(rows_run), rows_trial, string(rows_resp), ...
          'VariableNames', {'run','trial_idx','resp'});

if ~isempty(outCsv)
    writetable(T, outCsv);
    fprintf('Wrote table: %s\n', outCsv);
end

% --- helpers
function v = ifelse_first(vec, alt)
    if ~isempty(vec), v = vec(1); else, v = alt; end
end
end

function [lab, t] = i_read_label_time(inFile)
    try
        TT = readtable(inFile);
        if width(TT) >= 2
            lab = string(TT{:,1});
            t   = double(TT{:,2});
            return;
        end
    catch
    end
    fid = fopen(inFile,'r'); assert(fid>0, 'Cannot open %s', inFile);
    C = textscan(fid, '%s %f', 'Delimiter', {',',';','\t',' '}, ...
                 'MultipleDelimsAsOne', true, 'HeaderLines', 0);
    fclose(fid);
    lab = string(C{1});
    t   = double(C{2});
end
