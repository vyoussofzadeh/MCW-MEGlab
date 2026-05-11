function p = permutation_test_auc(X, Y, groups, evalFun, B, seed)
% evalFun: @(X,Y,groups) returns a single AUC on grouped CV
if nargin<5, B=200; end, if nargin<6, seed=0; end
rng(seed);
auc0 = evalFun(X,Y,groups); cnt=0;
for b=1:B
  Yperm = Y(randperm(numel(Y))); 
  aucb = evalFun(X,Yperm,groups);
  if aucb >= auc0, cnt=cnt+1; end
end
p = (cnt+1)/(B+1);
fprintf('Observed AUC=%.3f, permutation p=%.4f (B=%d)\n', auc0, p, B);
end
