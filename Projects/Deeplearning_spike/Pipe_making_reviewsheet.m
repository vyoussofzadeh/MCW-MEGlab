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

%%
%% === Assign overlapping reviewer pairs per SUBJECT (balanced + randomized) ===
% Build all unique reviewer pairs (for 4 reviewers -> 6 pairs)
pairs = nchoosek(1:numel(physicians), 2);   % [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]

% Subjects in stable order; we will assign pairs in randomized balanced blocks
[subs2, ~, g2] = unique(T.Subject, 'stable');
nSub = numel(subs2);

% --- Reproducible randomness (change seed if you want a new draw)
randSeed = 20250912; 
rng(randSeed);

% Build a randomized, balanced sequence of pairs of length nSub
pairSeq = zeros(nSub, 2);
writeIdx = 1;
while writeIdx <= nSub
    order = randperm(size(pairs,1));           % new random order of the 6 pairs
    take  = min(size(pairs,1), nSub - writeIdx + 1);
    pairSeq(writeIdx:writeIdx+take-1, :) = pairs(order(1:take), :);
    writeIdx = writeIdx + take;
end

% Optional: randomly flip P1/P2 per subject to avoid systematic labeling
flipMask = rand(nSub,1) > 0.5;
pairSeq(flipMask,:) = pairSeq(flipMask, [2 1]);

% Write assignments to all rows for each subject (both runs share same pair)
P1 = cell(height(T),1);
P2 = cell(height(T),1);
for s = 1:nSub
    idx = find(g2 == s);  % rows for this subject (2 runs)
    [P1{idx}] = deal(physicians{pairSeq(s,1)});
    [P2{idx}] = deal(physicians{pairSeq(s,2)});
end

T.Physician1 = P1;
T.Physician2 = P2;

% (Optional) quick load summary per reviewer (counts per *run*)
counts = struct(); for p = 1:numel(physicians), counts.(physicians{p}) = 0; end
for s = 1:nSub
    % +1 per run; 2 runs per subject -> +2 to each reviewer in the pair
    counts.(physicians{pairSeq(s,1)}) = counts.(physicians{pairSeq(s,1)}) + sum(g2==s);
    counts.(physicians{pairSeq(s,2)}) = counts.(physicians{pairSeq(s,2)}) + sum(g2==s);
end
disp(counts);

%% === Export ===
writetable(T, outCsv);
fprintf('Assigned %d subjects (2 runs each; %d files) -> %s\n', numel(subs2), height(T), outCsv);

