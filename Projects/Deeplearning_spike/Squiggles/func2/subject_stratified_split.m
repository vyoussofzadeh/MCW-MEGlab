function [trMaskSub, vaMaskSub] = subject_stratified_split(ySub, valFrac, seed)
% ySub: categorical labels for one subjects epochs (length Ns)
% Ensures both classes appear in TRAIN and VAL when counts allow.

if nargin<2 || isempty(valFrac), valFrac = 0.2; end
if nargin<3, seed = 0; end
rng(seed,'twister');

Ns = numel(ySub);
classes = categories(ySub);
trMaskSub = false(Ns,1);
vaMaskSub = false(Ns,1);

% Per-class indices
idxC = cell(numel(classes),1);
for c = 1:numel(classes)
    idxC{c} = find(ySub == classes{c});
end

% If only one class exists, do a simple holdout (cant stratify across classes)
onlyOneClass = sum(cellfun(@numel, idxC) > 0) == 1;
if onlyOneClass
    nVal = max(1, round(valFrac * Ns));
    perm = randperm(Ns);
    vaMaskSub(perm(1:nVal)) = true;
    trMaskSub(~vaMaskSub)   = true;
    return;
end

% Two-class (or more) stratified: hold out from each class
for c = 1:numel(classes)
    idx = idxC{c};
    nc  = numel(idx);
    if nc == 0, continue; end
    nValC = max(1, round(valFrac * nc));
    % ensure at least 1 remains in TRAIN if class has >=2
    if nc == 1
        % put that single sample into VAL only if other classes still have TRAIN+VAL;
        % otherwise keep it in TRAIN to avoid empty TRAIN.
        nValC = 0;
    elseif nc == 2
        nValC = 1;
    end
    sel = idx(randperm(nc, nValC));
    vaMaskSub(sel) = true;
end

trMaskSub = ~vaMaskSub;

% Safety: if any class missing in TRAIN or VAL, rebalance minimally
for side = 1:2
    mask = trMaskSub; name = 'TRAIN';
    if side==2, mask = vaMaskSub; name = 'VAL'; end
    for c = 1:numel(classes)
        if ~any(ySub(mask)==classes{c})
            % Move one sample of this class from the other side
            from = find(ySub(~mask)==classes{c});
            if ~isempty(from)
                if side==1
                    % need one in TRAIN -> move from VAL to TRAIN
                    idxFrom = find(vaMaskSub,1,'first');
                    vaMaskSub(idxFrom) = false; trMaskSub(idxFrom) = true;
                else
                    % need one in VAL -> move from TRAIN to VAL
                    idxFrom = find(trMaskSub & (ySub==classes{c}), 1, 'first');
                    if ~isempty(idxFrom)
                        trMaskSub(idxFrom) = false; vaMaskSub(idxFrom) = true;
                    end
                end
            end
        end
    end
end

% Final guards: ensure both sets non-empty
if ~any(trMaskSub), trMaskSub(1) = true; vaMaskSub(1) = false; end
if ~any(vaMaskSub)
    j = find(trMaskSub,1,'last');
    vaMaskSub(j) = true; trMaskSub(j) = false;
end
end
