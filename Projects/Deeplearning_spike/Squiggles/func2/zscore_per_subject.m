function [XtrZ, XvaZ] = zscore_per_subject(Xtr, Ytr, Gtr, Xva, Gva)
% Fit mean/std per subject on that subjects TRAIN epochs, apply to both his TRAIN & VAL
subs = unique(Gtr,'stable');
% build dict
muMap = containers.Map; sdMap = containers.Map;
for s = 1:numel(subs)
    ss = subs{s};
    rows = strcmp(Gtr, ss);
    if ~any(rows), continue; end
    % concatenate time across epochs (per channel row)
    xcat = cat(2, Xtr{rows});
    mu = mean(xcat, 2); sd = std(xcat, 0, 2); sd(sd<1e-6)=1;
    muMap(ss) = mu; sdMap(ss) = sd;
end
% apply
XtrZ = cell(size(Xtr)); XvaZ = cell(size(Xva));
for i=1:numel(Xtr)
    mu=muMap(Gtr{i}); sd=sdMap(Gtr{i}); XtrZ{i} = (Xtr{i} - mu) ./ sd;
end
for i=1:numel(Xva)
    ss = Gva{i};
    if isKey(muMap, ss)
        mu=muMap(ss); sd=sdMap(ss); XvaZ{i} = (Xva{i} - mu) ./ sd;
    else
        % unseen subject in VAL (rare if split by subjects) -> global z
        xcat = cat(2, Xtr{:}); mu=mean(xcat,2); sd=std(xcat,0,2); sd(sd<1e-6)=1;
        XvaZ{i} = (Xva{i} - mu) ./ sd;
    end
end
end
