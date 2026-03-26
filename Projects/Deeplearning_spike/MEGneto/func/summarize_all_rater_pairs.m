function [PairSummary, WorstRuns, FullRuns] = summarize_all_rater_pairs(labelsMat, base_cell, raters)
% labelsMat: N×R (0/1/NaN), base_cell: N×1 cellstr, raters: 1×R cellstr
% Outputs:
%   PairSummary: per reviewer-pair totals
%   WorstRuns : top disagreement bases per pair
%   FullRuns  : ALL bases per pair with overlap/disagreement/agreeRate/kappa

R = numel(raters);
pairs = nchoosek(1:R, 2);
bases = unique(base_cell);

% ---- overall pair summaries ----
out = cell(size(pairs,1), 11);
for p = 1:size(pairs,1)
    i = pairs(p,1); j = pairs(p,2);
    [nOverlap, nAgree1, nAgree0, nDis, agreeRate, kappa] = pair_stats(labelsMat(:,i), labelsMat(:,j));

    out(p,:) = {raters{i}, raters{j}, nOverlap, nAgree1, nAgree0, nDis, agreeRate, kappa, ...
                nDis / max(1,nOverlap), mean(labelsMat(~isnan(labelsMat(:,i)),i)), mean(labelsMat(~isnan(labelsMat(:,j)),j))};
end

PairSummary = cell2table(out, 'VariableNames', { ...
    'Rater1','Rater2','nOverlap','bothAgree','bothReject','disagree','agreeRate','kappa', ...
    'disagreeRate','pApprove_R1','pApprove_R2'});
PairSummary = sortrows(PairSummary, {'kappa','nOverlap'}, {'descend','descend'});

% ---- FULL per-run list for every pair ----
fullRows = {};
for p = 1:size(pairs,1)
    i = pairs(p,1); j = pairs(p,2);
    r1 = raters{i}; r2 = raters{j};

    for b = 1:numel(bases)
        idxB = strcmp(base_cell, bases{b});
        a = labelsMat(idxB,i);
        b2 = labelsMat(idxB,j);
        mask = ~isnan(a) & ~isnan(b2);
        nO = sum(mask);
        if nO == 0
            continue; % only store runs with overlap
        end

        aa = a(mask); bb = b2(mask);
        bothAgree  = sum(aa==1 & bb==1);
        bothReject = sum(aa==0 & bb==0);
        dis        = sum(aa~=bb);
        ar         = (bothAgree+bothReject)/nO;

        % per-run kappa
        pa1 = mean(aa==1); pb1 = mean(bb==1);
        pe  = pa1*pb1 + (1-pa1)*(1-pb1);
        kap = (ar - pe) / max(eps, (1 - pe));

        fullRows(end+1,:) = {r1, r2, bases{b}, nO, bothAgree, bothReject, dis, dis/nO, ar, kap}; %#ok<AGROW>
    end
end

FullRuns = cell2table(fullRows, 'VariableNames', ...
    {'Rater1','Rater2','Base','nOverlap','bothAgree','bothReject','disagreeN','disagreeRate','agreeRate','kappa'});

% sort so the most problematic are on top
FullRuns = sortrows(FullRuns, {'disagreeRate','nOverlap'}, {'descend','descend'});

% ---- WorstRuns: top 10 per pair (from FullRuns) ----
WorstRuns = {};
topN = 10;
minOverlapForListing = 5;

for p = 1:size(pairs,1)
    i = pairs(p,1); j = pairs(p,2);
    r1 = raters{i}; r2 = raters{j};

    idx = strcmp(FullRuns.Rater1, r1) & strcmp(FullRuns.Rater2, r2) & FullRuns.nOverlap >= minOverlapForListing;
    Tpair = FullRuns(idx,:);

    if isempty(Tpair), continue; end
    Tpair = sortrows(Tpair, {'disagreeRate','nOverlap'}, {'descend','descend'});
    Tpair = Tpair(1:min(topN,height(Tpair)), :);

    for t = 1:height(Tpair)
        WorstRuns(end+1,:) = {r1, r2, Tpair.Base{t}, Tpair.nOverlap(t), Tpair.disagreeN(t), Tpair.disagreeRate(t), Tpair.agreeRate(t)}; %#ok<AGROW>
    end
end

WorstRuns = cell2table(WorstRuns, 'VariableNames', ...
    {'Rater1','Rater2','Base','nOverlap','disagreeN','disagreeRate','agreeRate'});
WorstRuns = sortrows(WorstRuns, {'disagreeRate','nOverlap'}, {'descend','descend'});
end


function [nOverlap, bothAgree, bothReject, disagree, agreeRate, kappa] = pair_stats(a, b)
mask = ~isnan(a) & ~isnan(b);
aa = a(mask); bb = b(mask);

nOverlap  = numel(aa);
bothAgree = sum(aa==1 & bb==1);
bothReject= sum(aa==0 & bb==0);
disagree  = sum(aa~=bb);
agreeRate = (bothAgree+bothReject) / max(1,nOverlap);

if nOverlap == 0
    kappa = NaN; return;
end
po = agreeRate;
pa1 = mean(aa==1);
pb1 = mean(bb==1);
pe = pa1*pb1 + (1-pa1)*(1-pb1);
kappa = (po - pe) / max(eps, (1 - pe));
end
