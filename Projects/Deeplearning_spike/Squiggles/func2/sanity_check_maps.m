function sanity_check_maps(X2D, tag)
% X2D: cell of C×T×D
C = size(X2D{1},1); T = size(X2D{1},2); D = size(X2D{1},3);
bad = [];
for i=1:numel(X2D)
    s = size(X2D{i}); if numel(s)<3, s(3)=1; end
    if s(1)~=C || s(2)~=T || s(3)~=D || any(~isfinite(X2D{i}(:)))
        bad(end+1)=i; %#ok<AGROW>
    end
end
if isempty(bad)
    fprintf('[%s] OK: C=%d T=%d D=%d (N=%d)\n', tag, C,T,D,numel(X2D));
else
    fprintf('[%s] BAD %d/%d; first: %s\n', tag, numel(bad), numel(X2D), mat2str(bad(1:min(10,end))));
    error('Fix shapes/NaNs before training');
end
end