function Zcells = stacks_to_1Dfeatures(Xcells)
Zcells = cell(size(Xcells));
for i=1:numel(Xcells)
    X = Xcells{i};
    Zcells{i} = reshape(X, size(X,1)*size(X,3), size(X,2));  % (C*D)×T
end
end
