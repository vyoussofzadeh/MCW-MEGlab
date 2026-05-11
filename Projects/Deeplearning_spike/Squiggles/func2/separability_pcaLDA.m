function separability_pcaLDA(X, Y, groups, kFeat, seed)
% X: {N×1} C×T, Y: categorical('NoSpike','Spike'), groups: {N×1}
if nargin<4, kFeat=16; end, if nargin<5, seed=0; end
rng(seed,'twister');
C = size(X{1},1);

% Fit channel PCA on TRAIN folds only inside CV
cv = grouped_kfold(groups, 5, seed);   % 5-fold grouped CV (helper below)
AUC=[]; BA=[];

for f=1:numel(cv)
    tr=cv{f}.train; va=cv{f}.val;
    % --- channel PCA from TRAIN ---
    covm=zeros(C); n=0;
    for i=find(tr).'
        Xi=double(X{i}); covm = covm + (Xi*Xi.')/size(Xi,2); n=n+1;
    end
    covm=covm/n; [U,~]=svd(covm,'econ'); W=U(:,1:max(kFeat, C));   % many PCs now
    % epoch features = log-var of first r PCs (use smaller r=kFeat)
    feat = @(S,idx) cell2mat(cellfun(@(x) log(var((W.'*double(x)),0,2)+eps).', S(idx), 'uni',0));
    Ftr = feat(X, tr); Fva = feat(X, va);
    Ftr = Ftr(:,1:kFeat); Fva = Fva(:,1:kFeat);
    % z-score features using TRAIN
    mu=mean(Ftr,1); sd=std(Ftr,0,1)+eps; Ztr=(Ftr-mu)./sd; Zva=(Fva-mu)./sd;

    % --- LDA (linear) ---
    M = fitcdiscr(Ztr, Y(tr), 'DiscrimType','linear','Gamma',0.0,'FillCoeffs','off');
    Yhat = predict(M, Zva);

    % metrics
    [ba, auc] = metrics_BA_AUC(Zva, Y(va), M);
    BA(end+1)=ba; AUC(end+1)=auc; %#ok<AGROW>
end
fprintf('PCA(log-var) + LDA  |  CV BA=%.3f±%.3f  AUC=%.3f±%.3f\n', mean(BA),std(BA), mean(AUC),std(AUC));
end

function [ba, auc] = metrics_BA_AUC(Z, Y, model)
% Balanced accuracy from confusion; AUC from scores (LDA posterior for 'Spike')
order = categories(Y); cm = confusionmat(Y, predict(model,Z),'Order',order);
rec = diag(cm)./max(1,sum(cm,2)); ba = mean(rec);
if any(strcmp(order,'Spike')), [~,score] = resubPredict(model); %#ok<ASGLU>
    [~,scoreVa]= predict(model,Z); pos = strcmp(order,'Spike');
    sp = scoreVa(:,pos); [~,~,~,auc] = perfcurve(Y, sp, 'Spike'); else, auc=NaN; end
end

function cv = grouped_kfold(groups, K, seed)
% Returns cell array of K folds with .train/.val logical masks (group-wise)
rng(seed); u = unique(groups,'stable'); G=numel(u);
perm = u(randperm(G)); splits = round(linspace(0,G,K+1));
cv = cell(1,K);
for k=1:K
    vg = perm(splits(k)+1:splits(k+1));
    val = ismember(groups, vg); train = ~val;
    cv{k} = struct('train',train,'val',val);
end
end
