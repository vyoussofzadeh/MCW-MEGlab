function [Xc,Yc,Gc] = cap_by_group(X,Y,G,capPerClass)
uG = unique(G,'stable'); Xc={}; Yc=categorical; Gc={};
for gi = 1:numel(uG)
  idxG = find(G==uG{gi});
  for c = categories(Y).'
    idx = idxG(Y(idxG)==c);
    if isempty(idx), continue; end
    if numel(idx) > capPerClass
      idx = idx(randperm(numel(idx), capPerClass));
    end
    Xc = [Xc; X(idx)]; Yc = [Yc; Y(idx)]; Gc = [Gc; G(idx)]; %#ok<AGROW>
  end
end
[Yc] = reordercats(Yc, categories(Y));
end
