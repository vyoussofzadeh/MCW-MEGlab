function separability_CSP(X, Y, groups, Kpairs, seed)
if nargin<4, Kpairs=3; end, if nargin<5, seed=0; end
cv = grouped_kfold(groups,5,seed); AUC=[]; BA=[];
for f=1:numel(cv)
  tr=cv{f}.train; va=cv{f}.val; C=size(X{1},1);
  % composite covariances from TRAIN
  S0=zeros(C); S1=zeros(C); n0=0; n1=0;
  for i=find(tr).'
    Xi=double(X{i}); Xi=Xi-mean(Xi,2);
    S=(Xi*Xi.')/size(Xi,2);
    if Y(i)=='Spike', S1=S1+S; n1=n1+1; else, S0=S0+S; n0=n0+1; end
  end
  S0=S0/max(n0,1); S1=S1/max(n1,1); S=S0+S1;
  [E,D]=eig((S+S')/2); d=diag(D); Ww = diag(1./sqrt(max(d,eps))) * E';
  S0w=Ww*S0*Ww'; S1w=Ww*S1*Ww'; [U,~]=eig(S1w, S0w+eps*eye(C)); A=U'*Ww;
  idx=[1:Kpairs, C-Kpairs+1:C]; A=A(idx,:);           % 2K filters

  % epoch features = log-var of CSP components
  feat = @(S,idx) cell2mat(cellfun(@(x) log(var((A*double(x)),0,2)+eps).', S(idx), 'uni',0));
  Ftr=feat(X,tr); Fva=feat(X,va);                     % N×(2K)
  mu=mean(Ftr); sd=std(Ftr)+eps; Ztr=(Ftr-mu)./sd; Zva=(Fva-mu)./sd;

  M=fitcdiscr(Ztr, Y(tr), 'DiscrimType','linear','Gamma',0.0);
  [ba,auc]=metrics_BA_AUC(Zva,Y(va),M); BA(end+1)=ba; AUC(end+1)=auc; %#ok<AGROW>
end
fprintf('CSP(%d) + LDA       |  CV BA=%.3f±%.3f  AUC=%.3f±%.3f\n', Kpairs, mean(BA),std(BA), mean(AUC),std(AUC));
end
