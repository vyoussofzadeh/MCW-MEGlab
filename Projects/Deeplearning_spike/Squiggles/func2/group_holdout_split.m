function [trainMask,valMask] = group_holdout_split(groups, valFrac, seed)
% Splits by unique group IDs to prevent leakage between train/val
if nargin<3, seed = 0; end
rng(seed,'twister');

u = unique(groups,'stable');
G = numel(u);
if G == 0
    error('group_holdout_split:NoGroups','No groups provided.');
elseif G == 1
    % Put the single group entirely into TRAIN; keep VAL empty
    valMask = false(size(groups));
    trainMask = true(size(groups));
    return;
end

nVal = max(1, min(round(valFrac * G), G-1));
idx = randperm(G);
valGroups = u(idx(1:nVal));
valMask = ismember(groups, valGroups);
trainMask = ~valMask;
end
