function [Yg, Yhat, cats] = unify_cats(Yg, Yhat, catsRef)
% Make labels same categorical type and category order
if nargin < 3 || isempty(catsRef)
    catsRef = categories(categorical(Yg));
end
Yg   = categorical(Yg,   catsRef, 'Ordinal', false);
Yhat = categorical(Yhat, catsRef, 'Ordinal', false);
% Drop/merge any unknown predicted labels to first category
extra = setdiff(categories(Yhat), catsRef);
if ~isempty(extra)
    Yhat = mergecats(Yhat, extra, catsRef{1});
    Yhat = reordercats(Yhat, catsRef);
end
cats = catsRef;
end
