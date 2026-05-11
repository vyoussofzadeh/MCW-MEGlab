function Xpruned = prune_channels(Xcells, idxKeep)
% Apply the same channel subset to every epoch
Xpruned = cellfun(@(x) x(idxKeep, :), Xcells, 'UniformOutput', false);
end
