function [E, meta] = load_epochs_mat(fn)
%LOAD_EPOCHS_MAT Load epochs array from a spike/no-spike .mat
% Returns:
%   E    : epochs array (whatever orientation in file)
%   meta : struct with Fs/labels/time if present

S = load(fn);
meta = struct();

% Optional metadata (keeps what exists)
if isfield(S,'Fs');      meta.Fs = double(S.Fs); end
if isfield(S,'labels');  meta.labels = S.labels; end
if isfield(S,'channels');meta.channels = S.channels; end
if isfield(S,'time');    meta.time = S.time; end
if isfield(S,'t');       meta.t = S.t; end

% Common epoch variable names
candidates = {'SpikeEpochs','NotSpikeEpochs','NoSpikeEpochs','epochs','Epochs','X'};
E = [];

for i = 1:numel(candidates)
    nm = candidates{i};
    if isfield(S, nm) && isnumeric(S.(nm)) && ndims(S.(nm)) == 3
        E = S.(nm);
        meta.varname = nm;
        return;
    end
end

% Fallback: pick the largest 3D numeric array in the file
f = fieldnames(S);
bestSz = 0; bestNm = '';
for i = 1:numel(f)
    v = S.(f{i});
    if isnumeric(v) && ndims(v) == 3
        sz = numel(v);
        if sz > bestSz
            bestSz = sz;
            bestNm = f{i};
        end
    end
end

if ~isempty(bestNm)
    E = S.(bestNm);
    meta.varname = bestNm;
else
    error('No 3D numeric epochs array found in: %s', fn);
end
end
