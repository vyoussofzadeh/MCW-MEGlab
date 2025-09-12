function S = compute_global_indices_from_table(T, runOrder, outIdxTxt)
% compute_global_indices_from_table
%   From table T(run, trial_idx, resp), compute global 1..N indices for resp1/resp2
%   based on the specified runOrder (cellstr/strings of run base names as they appear in Brainstorm).
%
%   Returns struct S with:
%     .offsets (table of run, nTrials, offset)
%     .idx_resp1 (global indices)
%     .idx_resp2 (global indices)
%
%   If outIdxTxt is provided, writes a small text file with the two index lists.

if ischar(T) || isstring(T), T = readtable(T); end
T.run = string(T.run); T.resp = string(T.resp);

% Per-run counts
runs = string(runOrder(:));
nPerRun = zeros(numel(runs),1);
for i = 1:numel(runs)
    nPerRun(i) = sum(T.run==runs(i) & ismember(T.resp, ["resp1","resp2","timeout"]));
end
offsets = [0; cumsum(nPerRun(1:end-1))];
S.offsets = table(runs, nPerRun, offsets, 'VariableNames', {'run','nTrials','offset'});

% Global index for each row
globalIdx = zeros(height(T),1);
for i = 1:numel(runs)
    sel = T.run==runs(i);
    globalIdx(sel) = offsets(i) + T.trial_idx(sel);
end

S.idx_resp1 = sort(globalIdx(T.resp=="resp1"));
S.idx_resp2 = sort(globalIdx(T.resp=="resp2"));

if nargin>=3 && ~isempty(outIdxTxt)
    fid = fopen(outIdxTxt,'w'); assert(fid>0, 'Cannot write %s', outIdxTxt);
    fprintf(fid,'resp1:'); fprintf(fid,' %d', S.idx_resp1); fprintf(fid,'\n');
    fprintf(fid,'resp2:'); fprintf(fid,' %d', S.idx_resp2); fprintf(fid,'\n');
    fclose(fid);
end

% Print a quick summary
fprintf('Total trials (all runs): %d\n', sum(nPerRun));
fprintf('resp1: %d  |  resp2: %d  | timeouts: %d\n', ...
    sum(T.resp=="resp1"), sum(T.resp=="resp2"), sum(T.resp=="timeout"));
end
