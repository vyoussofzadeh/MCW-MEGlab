function [mdl, modelPCA, Ztr, Zva, cm, acc, bal] = train_pca_classifier(XTrain2, YTrain, XVal2, YVal, varargin)
% End-to-end: stats ? PCA ? logistic classifier ? metrics

p = inputParser;
addParameter(p,'fs',[],@(x) isempty(x) || (isscalar(x)&&x>0));
addParameter(p,'addBeta',false,@islogical);
addParameter(p,'var',95,@(x) isscalar(x) && x>0 && x<=100);
parse(p,varargin{:});
fs = p.Results.fs; addBeta = p.Results.addBeta; varexpl = p.Results.var;

% 1) Build per-epoch features
Ftr = build_epoch_stats(XTrain2, 'fs', fs, 'addBeta', addBeta);
Fva = build_epoch_stats(XVal2,   'fs', fs, 'addBeta', addBeta);

% 2) PCA (fit on TRAIN, apply to VAL)
[Ztr, modelPCA] = pca_fit_apply(Ftr, 'var', varexpl);
Zva = pca_apply(Fva, modelPCA);

% 3) Class weights
cats = categories(YTrain);
YTrain = reordercats(YTrain, cats);
YVal   = reordercats(YVal,   cats);
cnt = countcats(YTrain); w = sum(cnt)./max(cnt,1); w = w/mean(w);
Wtr = arrayfun(@(y) w(strcmp(cats, string(y))), YTrain);

% 4) Train a tiny classifier (logistic)
mdl = fitclinear(Ztr, YTrain, ...
    'Learner','logistic', 'Solver','lbfgs', 'ClassNames',cats, 'Weights',double(Wtr));

% 5) Evaluate
Yp = predict(mdl, Zva);
YValU = categorical(YVal, cats, 'Ordinal', false);
YpU   = categorical(Yp,   cats, 'Ordinal', false);
cm  = confusionmat(YValU, YpU, 'Order', cats);
acc = sum(diag(cm))/sum(cm,'all');
rec = diag(cm)./max(1,sum(cm,2));
bal = mean(rec);

fprintf('PCA+LogReg ? D=%d (%.1f%% var): Acc=%.3f  BA=%.3f\n', modelPCA.D, sum(modelPCA.expl(1:modelPCA.D)), acc, bal);
end
