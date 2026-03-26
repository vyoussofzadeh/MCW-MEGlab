function [RunPairStats, PairGroupRuns] = summarize_pairs_by_run(labelsMat, base_cell, raters)
% Summarize overlap/agreements at RUN level (Base), then aggregate per reviewer pair.
%
% Outputs:
%   RunPairStats  : one row per (Rater1,Rater2,Base)
%   PairGroupRuns : per pair micro + macro across runs

R = numel(raters);
pairs = nchoosek(1:R,2);
bases = unique(base_cell);

% -------- Run-level table ----------
rows = {};
for p = 1:size(pairs,1)
    i = pairs(p,1); j = pairs(p,2);
    r1 = raters{i}; r2 = raters{j};
    
    for b = 1:numel(bases)
        base = bases{b};
        idxB = strcmp(base_cell, base);
        
        a = labelsMat(idxB,i);
        c = labelsMat(idxB,j);
        
        %         [nO, bA, bR, dis, ar, kap, pA1, pB1] = pair_stats_with_rates(a, c);
        
        st = pair_stats_with_rates(a, c);
        
        if st.nOverlap > 0
            rows(end+1,:) = {r1,r2,base, ...
                st.nOverlap, st.bothAgree, st.bothReject, st.disagreeN, st.disagreeRate, st.agreeRate, st.kappa, ...
                st.r1_n1, st.r1_n0, st.r2_n1, st.r2_n0, ...
                st.pApprove_R1, st.pApprove_R2}; %#ok<AGROW>
        end
        
%         if nO > 0
%             rows(end+1,:) = {r1,r2,base,nO,bA,bR,dis,dis/max(1,nO),ar,kap,pA1,pB1}; %#ok<AGROW>
%         end
    end
end

% RunPairStats = cell2table(rows, 'VariableNames', { ...
%     'Rater1','Rater2','Base','nOverlap','bothAgree','bothReject','disagreeN', ...
%     'disagreeRate','agreeRate','kappa','pApprove_R1','pApprove_R2'});

RunPairStats = cell2table(rows, 'VariableNames', { ...
    'Rater1','Rater2','Base', ...
    'nOverlap','bothAgree','bothReject','disagreeN','disagreeRate','agreeRate','kappa', ...
    'R1_nAccept1','R1_nReject0','R2_nAccept1','R2_nReject0', ...
    'pApprove_R1','pApprove_R2'});

RunPairStats = sortrows(RunPairStats, {'Rater1','Rater2','nOverlap'}, {'ascend','ascend','descend'});

% -------- Group-level aggregation per pair (across runs) ----------
grows = {};
for p = 1:size(pairs,1)
    i = pairs(p,1); j = pairs(p,2);
    r1 = raters{i}; r2 = raters{j};
    
    Tp = RunPairStats(strcmp(RunPairStats.Rater1,r1) & strcmp(RunPairStats.Rater2,r2), :);
    if isempty(Tp), continue; end
    
    % MICRO: sum counts across runs then recompute
    nO = sum(Tp.nOverlap);
    bA = sum(Tp.bothAgree);
    bR = sum(Tp.bothReject);
    dis = sum(Tp.disagreeN);
    
    agreeRate_micro = (bA+bR)/max(1,nO);
    disagreeRate_micro = dis/max(1,nO);
    
    % pooled p(approve) approx weighted by overlap
    w = Tp.nOverlap / max(1,sum(Tp.nOverlap));
    pA = nansum(w .* Tp.pApprove_R1);
    pB = nansum(w .* Tp.pApprove_R2);
    pe = pA*pB + (1-pA)*(1-pB);
    kappa_micro = (agreeRate_micro - pe) / max(eps, (1-pe));
    
    % MACRO: average run-wise stats equally (each run = 1 vote)
    agreeRate_macro = mean(Tp.agreeRate, 'omitnan');
    disagreeRate_macro = mean(Tp.disagreeRate, 'omitnan');
    kappa_macro = mean(Tp.kappa, 'omitnan');
    
    nRunsOverlap = height(Tp);
    
    grows(end+1,:) = {r1,r2,nRunsOverlap,nO, ...
        agreeRate_micro, disagreeRate_micro, kappa_micro, ...
        agreeRate_macro, disagreeRate_macro, kappa_macro}; %#ok<AGROW>
end

PairGroupRuns = cell2table(grows, 'VariableNames', { ...
    'Rater1','Rater2','nRunsWithOverlap','nOverlapCandidates', ...
    'agreeRate_micro','disagreeRate_micro','kappa_micro', ...
    'agreeRate_macro','disagreeRate_macro','kappa_macro'});

PairGroupRuns = sortrows(PairGroupRuns, {'kappa_micro','nOverlapCandidates'}, {'descend','descend'});
end

function stats = pair_stats_with_rates(a, b)
% a,b: vectors with {0,1,NaN} for ONE run and ONE rater pair
% Returns a struct with overlap-only counts and rates.

mask = ~isnan(a) & ~isnan(b);
aa = a(mask); bb = b(mask);

stats.nOverlap   = numel(aa);

stats.r1_n1 = sum(aa==1);
stats.r1_n0 = sum(aa==0);
stats.r2_n1 = sum(bb==1);
stats.r2_n0 = sum(bb==0);

stats.bothAgree  = sum(aa==1 & bb==1);
stats.bothReject = sum(aa==0 & bb==0);
stats.disagreeN  = sum(aa~=bb);

stats.agreeRate    = (stats.bothAgree + stats.bothReject) / max(1, stats.nOverlap);
stats.disagreeRate = stats.disagreeN / max(1, stats.nOverlap);

if stats.nOverlap == 0
    stats.pApprove_R1 = NaN; stats.pApprove_R2 = NaN; stats.kappa = NaN;
    return;
end

stats.pApprove_R1 = stats.r1_n1 / stats.nOverlap;
stats.pApprove_R2 = stats.r2_n1 / stats.nOverlap;

pe = stats.pApprove_R1*stats.pApprove_R2 + (1-stats.pApprove_R1)*(1-stats.pApprove_R2);
stats.kappa = (stats.agreeRate - pe) / max(eps, (1 - pe));
end

% function [nOverlap, bothAgree, bothReject, disagree, agreeRate, kappa, pA1, pB1] = pair_stats_with_rates(a, b)
% mask = ~isnan(a) & ~isnan(b);
% aa = a(mask); bb = b(mask);
%
% nOverlap  = numel(aa);
% bothAgree = sum(aa==1 & bb==1);
% bothReject= sum(aa==0 & bb==0);
% disagree  = sum(aa~=bb);
% agreeRate = (bothAgree+bothReject) / max(1,nOverlap);
%
% if nOverlap == 0
%     kappa = NaN; pA1 = NaN; pB1 = NaN; return;
% end
%
% pA1 = mean(aa==1);
% pB1 = mean(bb==1);
%
% pe = pA1*pB1 + (1-pA1)*(1-pB1);
% kappa = (agreeRate - pe) / max(eps, (1 - pe));
% end
