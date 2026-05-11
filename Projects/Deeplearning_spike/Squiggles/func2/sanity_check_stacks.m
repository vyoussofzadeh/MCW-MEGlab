function sanity_check_stacks(Xcells, tag)
C = size(Xcells{1},1); T = size(Xcells{1},2); D = size(Xcells{1},3);
bad = [];
for i=1:numel(Xcells)
    s = size(Xcells{i}); if numel(s)<3, s(3)=1; end
    if s(1)~=C || s(2)~=T || s(3)~=D || ~isfinite(sum(Xcells{i}(:)))
        bad(end+1) = i; %#ok<AGROW>
    end
end
if isempty(bad)
    fprintf('[%s] OK: C=%d, T=%d, D=%d across %d epochs\n', tag,C,T,D,numel(Xcells));
else
    fprintf('[%s] BAD %d/%d; first 10: %s\n', tag, numel(bad), numel(Xcells), mat2str(bad(1:min(10,end))));
    error('Fix shapes/NaNs before training');
end
end
