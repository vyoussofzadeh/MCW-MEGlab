function [SubjPairRuns, PairGroup] = summarize_pairs_by_subject(labelsMat, base_cell, raters)
% Summarize overlap/agreements at SUBJECT level, then aggregate per reviewer pair.
%
% Inputs:
%   labelsMat : N×R (0/1/NaN)
%   base_cell : N×1 cellstr, each entry is a "Base" string
%   raters    : 1×R cellstr
%
% Outputs:
%   SubjPairRuns : table with one row per (Rater1,Rater2,Subject)
%   PairGroup    : table with micro/macro summaries per pair

R = numel(raters);
pairs = nchoosek(1:R,2);

% Build subject_id per candidate from base_cell (strip _Run\d+ suffix)
subj_id = cellfun(@extract_subject_id, base_cell, 'UniformOutput', false);
subjects = unique(subj_id);

% -------- Subject-level table ----------
rows = {};
for p = 1:size(pairs,1)
    i = pairs(p,1); j = pairs(p,2);
    r1 = raters{i}; r2 = raters{j};

    for s = 1:numel(subjects)
        sid = subjects{s};
        idxS = strcmp(subj_id, sid);

        a = labelsMat(idxS,i);
        b = labelsMat(idxS,j);

        [nO, bA, bR, dis, ar, kap, pA1, pB1] = pair_stats_with_rates(a,b);

        if nO > 0
            rows(end+1,:) = {r1,r2,sid,nO,bA,bR,dis,dis/max(1,nO),ar,kap,pA1,pB1}; %#ok<AGROW>
        end
    end
end

SubjPairRuns = cell2table(rows, 'VariableNames', { ...
    'Rater1','Rater2','Subject','nOverlap','bothAgree','bothReject','disagreeN', ...
    'disagreeRate','agreeRate','kappa','pApprove_R1','pApprove_R2'});

% Sort for readability
SubjPairRuns = sortrows(SubjPairRuns, {'Rater1','Rater2','nOverlap'}, {'ascend','ascend','descend'});

% -------- Group-level aggregation per pair ----------
grows = {};
for p = 1:size(pairs,1)
    i = pairs(p,1); j = pairs(p,2);
    r1 = raters{i}; r2 = raters{j};

    Tp = SubjPairRuns(strcmp(SubjPairRuns.Rater1,r1) & strcmp(SubjPairRuns.Rater2,r2), :);
    if isempty(Tp)
        continue;
    end

    % MICRO: sum counts across subjects then recompute
    nO = sum(Tp.nOverlap);
    bA = sum(Tp.bothAgree);
    bR = sum(Tp.bothReject);
    dis = sum(Tp.disagreeN);
    agreeRate_micro = (bA+bR)/max(1,nO);
    disagreeRate_micro = dis/max(1,nO);

    % micro kappa: recompute from pooled label rates approx
    % We'll compute pooled p(approve) for each rater on overlapped items using subject-level rates weighted by nOverlap.
    % This is a good approximation and works well in practice.
    w = Tp.nOverlap / max(1,sum(Tp.nOverlap));
    pA = nansum(w .* Tp.pApprove_R1);
    pB = nansum(w .* Tp.pApprove_R2);
    pe = pA*pB + (1-pA)*(1-pB);
    kappa_micro = (agreeRate_micro - pe) / max(eps, (1-pe));

    % MACRO: average subject-wise stats (each subject equal weight)
    agreeRate_macro = mean(Tp.agreeRate, 'omitnan');
    disagreeRate_macro = mean(Tp.disagreeRate, 'omitnan');
    kappa_macro = mean(Tp.kappa, 'omitnan');

    nSubjectsOverlap = height(Tp);

    grows(end+1,:) = {r1,r2,nSubjectsOverlap,nO, ...
        agreeRate_micro, disagreeRate_micro, kappa_micro, ...
        agreeRate_macro, disagreeRate_macro, kappa_macro}; %#ok<AGROW>
end

PairGroup = cell2table(grows, 'VariableNames', { ...
    'Rater1','Rater2','nSubjectsWithOverlap','nOverlapCandidates', ...
    'agreeRate_micro','disagreeRate_micro','kappa_micro', ...
    'agreeRate_macro','disagreeRate_macro','kappa_macro'});

PairGroup = sortrows(PairGroup, {'kappa_micro','nOverlapCandidates'}, {'descend','descend'});
end


function [nOverlap, bothAgree, bothReject, disagree, agreeRate, kappa, pA1, pB1] = pair_stats_with_rates(a, b)
% a,b: vectors with {0,1,NaN} for ONE subject and ONE rater pair
mask = ~isnan(a) & ~isnan(b);
aa = a(mask); bb = b(mask);

nOverlap  = numel(aa);
bothAgree = sum(aa==1 & bb==1);
bothReject= sum(aa==0 & bb==0);
disagree  = sum(aa~=bb);
agreeRate = (bothAgree+bothReject) / max(1,nOverlap);

if nOverlap == 0
    kappa = NaN; pA1 = NaN; pB1 = NaN; return;
end

pA1 = mean(aa==1);
pB1 = mean(bb==1);

pe = pA1*pB1 + (1-pA1)*(1-pB1);
kappa = (agreeRate - pe) / max(eps, (1 - pe));
end


function sid = extract_subject_id(base)
% base: e.g. 'bednar_peggy_Run02_spont_supine_raw_t_sss_ecgClean'
tok = regexp(base, '^(.*)_Run\d+', 'tokens', 'once');
if ~isempty(tok), sid = tok{1};
else, sid = base;
end
end
