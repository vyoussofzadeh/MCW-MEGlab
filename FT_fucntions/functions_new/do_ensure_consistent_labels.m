function anot_data = do_ensure_consistent_labels(anot_data, refLabels)
% ENSURE_CONSISTENT_LABELS  Harmonise anot_data.label with trial matrix size
%
%   anot_data = ensure_consistent_labels(anot_data)
%       - trims   or pads channels so that numel(label)==size(trial,1)
%
%   anot_data = ensure_consistent_labels(anot_data, refLabels)
%       - additionally re-orders channels to exactly match refLabels
%         (missing channels are dropped; extra channels are ignored)
%
% Author: MCW MEG group  2025-05-02

% ---------- sanity check ----------
if ~iscell(anot_data.label)
    error('anot_data.label must be a cell array of channel names.');
end
if isempty(anot_data.trial) || ~iscell(anot_data.trial)
    error('anot_data.trial must be a non-empty cell array.');
end

% ---------- basic size alignment ----------
nLabel = numel(anot_data.label);
nChan  = size(anot_data.trial{1},1);

% Case 1: mismatch  keep common subset (min)
if nLabel ~= nChan
    warning(['Label/channel mismatch (labels=%d, chans=%d)  ' ...
             'keeping common subset.'], nLabel, nChan);
    keepIdx          = 1:min(nLabel,nChan);
    anot_data.label  = anot_data.label(keepIdx);
    for k = 1:numel(anot_data.trial)
        anot_data.trial{k} = anot_data.trial{k}(keepIdx,:);
    end
end

% ---------- optional re-order to reference ----------
if nargin > 1 && ~isempty(refLabels)
    [found, idxData] = ismember(refLabels, anot_data.label);

    % Drop channels absent from this epoch
    idxData(~found) = [];     %#ok<NASGU>
    refLabels(~found) = [];

    % Re-order
    anot_data.label  = refLabels;
    for k = 1:numel(anot_data.trial)
        anot_data.trial{k} = anot_data.trial{k}( idxData(found), : );
    end
end

% final assertion
assert(numel(anot_data.label) == size(anot_data.trial{1},1), ...
       'Labels still inconsistent after fix.');

end
