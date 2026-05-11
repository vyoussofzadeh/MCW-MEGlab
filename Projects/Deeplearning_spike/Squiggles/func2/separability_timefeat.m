function separability_timefeat(X, Y, groups, kPCA, seed)
if nargin<4, kPCA=32; end, if nargin<5, seed=0; end
cv = grouped_kfold(groups,5,seed); AUC=[]; BA=[];
featFun = @(x) [sum(abs(diff(x,1,2)),2); var(x,0,2); kurtosis(x,0,2)]; % 3C×1
for f=1:numel(cv)
  tr=cv{f}.train; va=cv{f}.val;
  Ftr = cell2mat(cellfun(@(x) featFun(double(x)), X(tr), 'uni',0)).'; % N×(3C)
  Fva = cell2mat(cellfun(@(x) featFun(double(x)), X(va), 'uni',0)).';
  [Cmat, Ztr, ~] = pca( (Ftr - mean(Ftr,1))./ (std(Ftr,0,1)+eps), 'NumComponents', kPCA);
  Zva = (Fva - mean(Ftr,1))./(std(Ftr,0,1)+eps) * Cmat;
  M = fitcsvm(Ztr, Y(tr), 'KernelFunction','linear', 'Standardize',true, ...
              'ClassNames',categories(Y), 'BoxConstraint',1);
  [ba,auc]=metrics_BA_AUC(Zva,Y(va),MdlToLDA(M,Ztr,Y(tr))); BA(end+1)=ba; AUC(end+1)=auc; %#ok<AGROW>
end
fprintf('Time feats + linSVM |  CV BA=%.3f±%.3f  AUC=%.3f±%.3f\n', mean(BA),std(BA), mean(AUC),std(AUC));
end

function L = MdlToLDA(SVM, Z, Y)
% Convert linear SVM to a linear score model for AUC/BA convenience
w = SVM.Beta; b = SVM.Bias;  % decision w'*x + b
L = fitcdiscr(Z, Y, 'DiscrimType','pseudoLinear'); %#ok<NASGU>
L.Coeffs(1,2).Linear = w; L.Coeffs(1,2).Const = b;  % rough plug-in
end
