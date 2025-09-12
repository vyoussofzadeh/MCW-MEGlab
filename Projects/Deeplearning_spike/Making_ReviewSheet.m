%%
% === Config ===
dirPath    = '/MEG_data/AHW_SpikeAnalysis/MEG_data/sss/Spike_times';
physicians = {'Manoj','Adi','Pradeep','Josh'};  % <- replace with real names
% outCsv     = fullfile(dirPath, 'review_assignments_first2runs.csv');

Savepath = '/MEG_data/AHW_SpikeAnalysis/Review_sheet';
outCsv     = fullfile(Savepath, 'review_assignments.csv');


%%
% Selection rule for which two runs per subject to keep:
% 'lowest'  -> two smallest run numbers (default)
% 'highest' -> two largest run numbers
% 'random'  -> two random runs (reproducible)
pickRule   = 'lowest';

%% === Collect all .eve files and parse Subject/Run ===
L = dir(fullfile(dirPath, '*.eve'));
S = struct('Subject',{},'RunNum',{},'RunLabel',{},'FileName',{},'FullPath',{});

for i = 1:numel(L)
    fname = L(i).name;
    % Parse "subject_runXX_" (case-insensitive, accepts Run6/run06/etc.)
    m = regexpi(fname, '^(?<subject>.+?)_run?0?(?<run>\d+)_', 'names', 'once');
    if isempty(m), continue; end
    rnum = str2double(m.run);
    if isnan(rnum), continue; end

    S(end+1) = struct( ... %#ok<SAGROW>
        'Subject', m.subject, ...
        'RunNum',  rnum, ...
        'RunLabel', sprintf('Run%02d', rnum), ...
        'FileName', fname, ...
        'FullPath', fullfile(dirPath, fname));
end

if isempty(S)
    error('No .eve files parsed in: %s', dirPath);
end

T = struct2table(S);

%% === Keep exactly TWO runs per subject (require subjects with >= 2 runs) ===
% Sort once to make selection deterministic
T = sortrows(T, {'Subject','RunNum'});

subs = unique(T.Subject, 'stable');
keepIdx = false(height(T),1);

switch lower(pickRule)
    case 'lowest'
        for s = 1:numel(subs)
            idx = find(strcmp(T.Subject, subs{s}));
            if numel(idx) < 2, continue; end
            % already sorted ascending by RunNum from sortrows
            keepIdx(idx(1:2)) = true;
        end
    case 'highest'
        for s = 1:numel(subs)
            idx = find(strcmp(T.Subject, subs{s}));
            if numel(idx) < 2, continue; end
            keepIdx(idx(end-1:end)) = true;
        end
    case 'random'
        rng(42); % reproducible selection
        for s = 1:numel(subs)
            idx = find(strcmp(T.Subject, subs{s}));
            if numel(idx) < 2, continue; end
            pick = randsample(idx, 2);
            keepIdx(pick) = true;
        end
    otherwise
        error('Unknown pickRule: %s', pickRule);
end

T = T(keepIdx,:);
if isempty(T)
    error('All rows were filtered out. Check parsing and that subjects have >= 2 runs.');
end

%% === Assign overlapping reviewer pairs per SUBJECT (both runs share the same pair) ===
% Pair cycle: [A,B] -> [B,C] -> [C,D] -> [D,A] -> repeat
pairCycle = [1 2; 2 3; 3 4; 4 1];

[subs2, ~, g2] = unique(T.Subject, 'stable');
P1 = cell(height(T),1);
P2 = cell(height(T),1);

for s = 1:numel(subs2)
    cyc = mod(s-1, size(pairCycle,1)) + 1;
    r1  = physicians{pairCycle(cyc,1)};
    r2  = physicians{pairCycle(cyc,2)};
    idx = find(g2 == s);           % both runs get same pair
    [P1{idx}] = deal(r1);
    [P2{idx}] = deal(r2);
end

T.Physician1 = P1;
T.Physician2 = P2;

% Nice ordering in the CSV
T = sortrows(T, {'Subject','RunNum'});

%% === Export ===
writetable(T, outCsv);
fprintf('Assigned %d subjects (2 runs each; %d files) -> %s\n', numel(subs2), height(T), outCsv);

