function cv = grouped_kfold(groups, K, seed)
% GROUPED_KFOLD  K-fold CV with group integrity (no subject/session leakage).
%   cv = grouped_kfold(groups, K, seed)
%   groups : {N×1} cell array of group IDs (e.g., subject/session)
%   K      : number of folds (e.g., 5)
%   seed   : RNG seed for reproducibility
% Returns:
%   cv : 1×K cell; each cv{k} has logical vectors .train and .val (length N)

if nargin<2 || isempty(K),    K = 5; end
if nargin<3 || isempty(seed), seed = 0; end
rng(seed,'twister');

N = numel(groups);
assert(iscell(groups) && N>0, 'groups must be {N×1} cell.');

[uG,~,gi] = unique(groups,'stable');   % gi maps sample -> group index
G = numel(uG);
K = max(1, min(K, G));                 % cannot have more folds than groups

% randomize group order, then slice into K folds as evenly as possible
permG = uG(randperm(G));
edges = round(linspace(0, G, K+1));
cv = cell(1,K);

for k = 1:K
    valGroups = permG(edges(k)+1 : edges(k+1));
    valMask   = ismember(groups, valGroups);
    trainMask = ~valMask;

    % if a fold accidentally becomes empty (can happen with tiny G), fix it
    if ~any(valMask)
        % move one random group from train to val
        trGroups = uG(~ismember(uG, valGroups));
        valGroups = trGroups(randi(numel(trGroups)));
        valMask   = ismember(groups, valGroups);
        trainMask = ~valMask;
    end

    cv{k} = struct('train', trainMask, 'val', valMask, ...
                   'valGroupIDs', {valGroups});
end
end
