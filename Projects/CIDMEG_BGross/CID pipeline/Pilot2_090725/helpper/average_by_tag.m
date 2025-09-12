function avgFile = average_by_tag(runNum, tagStr, outName)
% Average all trials in the current folder that match a given Run and tag in Comment.
% Usage: average_by_tag(2,'RT=fast','avg_run2_RTfast.mat')

% files = dir(sprintf('data_%d_trial*.mat', runNum));
files = dir('data_2_trial*.mat');
Fs = []; lastS = [];
keep = [];

for k = 1:numel(files)
    S = load(files(k).name);
    if isfield(S,'Comment') && contains(S.Comment, tagStr)
        Fs = cat(3, Fs, S.F);
        keep(end+1) = k; %#ok<AGROW>
        lastS = S;
    end
end

if isempty(keep), error('No files matched Run=%d and tag "%s".', runNum, tagStr); end

Favg = mean(Fs, 3);
Savg = lastS;
Savg.F       = Favg;
Savg.nAvg    = size(Fs,3);
Savg.Comment = sprintf('Avg |Run=%d| |%s| (N=%d)', runNum, tagStr, Savg.nAvg);

if isfield(Savg,'History') && iscell(Savg.History)
    Savg.History(end+1,1:3) = {datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
                               'average_by_tag', sprintf('Run=%d, Tag=%s, N=%d', runNum, tagStr, Savg.nAvg)};
end

if nargin < 3 || isempty(outName)
    outName = sprintf('data_%d_avg_%s.mat', runNum, regexprep(tagStr,'[^\w]+',''));
end
save(outName, '-struct', 'Savg');
fprintf('Wrote %s\n', outName);
avgFile = outName;
end
