clear

%% -------------------- META-ONLY: build decisions + summaries -------------
root   = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files';
raters = {'Adi','Josh','Manoj','Pradeep'};
Fs_default = 200;

% If sample indices might be slightly different across reviewers, set tolSec>0
% tolSec = 0      -> exact match on sample index
% tolSec = 0.010  -> 10 ms binning (requires Fs)
tolSec = 0;   % start with exact matching

spath = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/MEGneto/metaprocessed';
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/MEGneto/func')

%% 1) Load ALL meta decisions into one table: (Rater, Base, sample, code, agree)
All = table();
for r = 1:numel(raters)
    rname = raters{r};
    metaDir = fullfile(root, rname, 'meta_files');
    ff = dir(fullfile(metaDir, ['*_' rname '_meta_data.txt']));

    for k = 1:numel(ff)
        meta_fn = fullfile(ff(k).folder, ff(k).name);
        base = erase(ff(k).name, ['_' rname '_meta_data.txt']);

        meta = read_meta_txt(meta_fn, Fs_default, false);
        Tm = meta.table;
        if isempty(Tm) || ~all(ismember({'sample','agree','code'}, Tm.Properties.VariableNames))
            fprintf('Skipping meta (bad format): %s\n', meta_fn);
            continue
        end

        Tr = table();
        Tr.Rater = repmat({rname}, height(Tm), 1);
        Tr.Base  = repmat({base},  height(Tm), 1);
        Tr.Subject = repmat({extract_subject_id(base)}, height(Tm), 1);

        Tr.sample = round(Tm.sample);
        Tr.code   = round(Tm.code);
        Tr.agree  = round(Tm.agree);   % 1=approved spike, 0=rejected

        All = [All; Tr]; %#ok<AGROW>
    end
end

fprintf('Meta-only: total rows loaded = %d\n', height(All));
fprintf('Meta-only: unique bases (runs) = %d\n', numel(unique(All.Base)));

%% 2) Build keys and pivot to Y: [nKeys x nRaters] with NaN when not reviewed
% Optional tolerance via sample binning
Fs = Fs_default;
if tolSec > 0
    tolSamp = max(1, round(tolSec * Fs));
    All.sample_key = tolSamp * round(double(All.sample)/tolSamp);
else
    All.sample_key = All.sample;
end

All.key = strcat(All.Base, "|", string(All.sample_key), "|", string(All.code));

uKeys = unique(All.key);
R = numel(raters);

Y = NaN(numel(uKeys), R);   % labels per rater
baseCol = cell(numel(uKeys),1);
subjCol = cell(numel(uKeys),1);
sampCol = zeros(numel(uKeys),1);
codeCol = zeros(numel(uKeys),1);

for i = 1:numel(uKeys)
    idx = (All.key == uKeys(i));
    baseCol{i} = All.Base{find(idx,1,'first')};
    subjCol{i} = All.Subject{find(idx,1,'first')};
    sampCol(i) = All.sample_key(find(idx,1,'first'));
    codeCol(i) = All.code(find(idx,1,'first'));

    for rr = 1:R
        rname = raters{rr};
        jdx = idx & strcmp(All.Rater, rname);
        if any(jdx)
            % if duplicates exist, take the first (should not happen ideally)
            Y(i,rr) = All.agree(find(jdx,1,'first'));
        end
    end
end

%% 3) Meta-only coverage counts per rater
for r = 1:R
    has = ~isnan(Y(:,r));
    fprintf('%s: candidates reviewed=%d, unique bases=%d\n', ...
        raters{r}, sum(has), numel(unique(baseCol(has))));
end

%% 4) Run-level + group-level pair summaries (meta-only)
[RunPairStats_meta, PairGroupRuns_meta] = summarize_pairs_by_run(Y, baseCol, raters);

disp('=== META-only RunPairStats (pair x base/run) ===');
disp(RunPairStats_meta(1:min(30,height(RunPairStats_meta)),:));

disp('=== META-only PairGroup (micro+macro across runs) ===');
disp(PairGroupRuns_meta);

writetable(RunPairStats_meta,  fullfile(spath,'META_PairRunLevel.csv'));
writetable(PairGroupRuns_meta, fullfile(spath,'META_PairGroupLevel_byRun.csv'));

%% 5) Build meta-only consensus labels per candidate key
minRatersForConsensus = 2;
consensusRule = 'majority'; % ties -> NaN

nReviewed = sum(~isnan(Y),2);
keep = (nReviewed >= minRatersForConsensus);

cons = NaN(size(Y,1),1);
agreeCount = NaN(size(Y,1),1);
rejectCount = NaN(size(Y,1),1);

for i = 1:size(Y,1)
    labs = Y(i, ~isnan(Y(i,:)));
    if isempty(labs), continue; end
    agreeCount(i)  = sum(labs==1);
    rejectCount(i) = sum(labs==0);

    switch lower(consensusRule)
        case 'majority'
            m = mean(labs);
            if m > 0.5, cons(i)=1;
            elseif m < 0.5, cons(i)=0;
            else, cons(i)=NaN; % tie
            end
    end
end

keep = keep & ~isnan(cons);

Consensus = table();
Consensus.Base    = baseCol(keep);
Consensus.Subject = subjCol(keep);
Consensus.sample  = sampCol(keep);
Consensus.code    = codeCol(keep);
Consensus.nReviewed = nReviewed(keep);
Consensus.agreeCount = agreeCount(keep);
Consensus.rejectCount = rejectCount(keep);
Consensus.consensus = cons(keep);

% Save full rater matrix too (useful)
for rr = 1:R
    Consensus.(sprintf('agree_%s', raters{rr})) = Y(keep,rr);
end

fprintf('Consensus candidates kept = %d\n', height(Consensus));
writetable(Consensus, fullfile(spath, 'META_ConsensusCandidates.csv'));

% Optionally save as MAT for downstream
save( fullfile(spath,'META_only_decisions.mat'), 'All', 'Y', 'baseCol', 'subjCol', 'sampCol', 'codeCol', 'Consensus', 'raters', ...
     'minRatersForConsensus', 'consensusRule', 'tolSec', '-v7.3');

%% 6) Rate reviewers from META-only decisions (reliability + coverage)
ReviewerRating = rate_reviewers_from_meta(Y, baseCol, raters);

disp('=== Reviewer rating (meta-only) ===');
disp(ReviewerRating);

writetable(ReviewerRating, fullfile(spath,'META_ReviewerRating.csv'));
 
 
%% -------- helper: subject id from base --------
function sid = extract_subject_id(base)
    tok = regexp(base, '^(.*)_Run\d+', 'tokens', 'once');
    if ~isempty(tok), sid = tok{1}; else, sid = base; end
end


function ReviewerRating = rate_reviewers_from_meta(Y, baseCol, raters)
% Rate reviewers using meta-only labels matrix Y (nKeys x nRaters).
% Outputs per-rater:
%   - coverage (#candidates, #bases)
%   - mean weighted kappa vs others
%   - mean weighted agreement rate vs others

R = numel(raters);

% Coverage per rater
nReviewed = sum(~isnan(Y),1)';  % R x 1
nBases = zeros(R,1);
for rr=1:R
    has = ~isnan(Y(:,rr));
    nBases(rr) = numel(unique(baseCol(has)));
end

% Pairwise kappas + agreement + overlaps
K = NaN(R,R);
A = NaN(R,R);
O = zeros(R,R);

for i=1:R
    for j=i+1:R
        a = Y(:,i); b = Y(:,j);
        mask = ~isnan(a) & ~isnan(b);
        aa = a(mask); bb = b(mask);
        nO = numel(aa);
        O(i,j) = nO; O(j,i) = nO;

        if nO == 0
            continue;
        end

        agreeRate = mean(aa==bb);
        pa1 = mean(aa==1); pb1 = mean(bb==1);
        pe = pa1*pb1 + (1-pa1)*(1-pb1);
        kappa = (agreeRate - pe) / max(eps, (1-pe));

        K(i,j) = kappa; K(j,i) = kappa;
        A(i,j) = agreeRate; A(j,i) = agreeRate;
    end
end

% Per-rater weighted averages vs others
meanK = NaN(R,1);
meanA = NaN(R,1);
totOverlap = zeros(R,1);

for i=1:R
    w = O(i,:);           % overlaps with others
    k = K(i,:);
    a = A(i,:);

    valid = (w>0) & ~isnan(k);
    if any(valid)
        meanK(i) = sum(w(valid).*k(valid)) / sum(w(valid));
        meanA(i) = sum(w(valid).*a(valid)) / sum(w(valid));
        totOverlap(i) = sum(w(valid));
    else
        meanK(i) = NaN;
        meanA(i) = NaN;
        totOverlap(i) = 0;
    end
end

% Optional: a single composite score
% (weights reliability more than coverage; tune as you like)
covScore = log10(1 + nReviewed);                % diminishing returns
relScore = max(0, meanK);                       % don't reward negative kappas
score = relScore + 0.15 * covScore;             % composite

ReviewerRating = table(raters(:), nReviewed, nBases, totOverlap, meanA, meanK, score, ...
    'VariableNames', {'Rater','nCandidatesReviewed','nBasesReviewed','totalPairOverlap', ...
                      'meanAgreeRate_weighted','meanKappa_weighted','score'});

ReviewerRating = sortrows(ReviewerRating, 'score', 'descend');
end
