function [Xtr_tem, Xva_tem, inputSize, Tseq] = timelock_to_seq1d(tl_tr, tl_va, varargin)
% Convert FieldTrip timelock (rpt×chan×time) to cell arrays for 1D models.
% Name-Value:
%   'UsePCA'  (false/true)  apply PCA across channels using TRAIN only
%   'VarKeep' (98..99)      % variance to keep if PCA is used
p = inputParser;
p.addParameter('UsePCA', false, @(x)islogical(x)||isscalar(x));
p.addParameter('VarKeep', 98, @isscalar);
p.parse(varargin{:});
usePCA = p.Results.UsePCA; varKeep = p.Results.VarKeep;

Ntr = size(tl_tr.trial,1); C = size(tl_tr.trial,2); T = size(tl_tr.trial,3);
Nva = size(tl_va.trial,1);

% --- pack to cell arrays [C × T]
Xtr = cell(Ntr,1);  for i=1:Ntr, Xtr{i} = squeeze(tl_tr.trial(i,:,:)); end
Xva = cell(Nva,1);  for i=1:Nva, Xva{i} = squeeze(tl_va.trial(i,:,:)); end

if usePCA
    % Fit PCA on TRAIN across channels; concatenate time across trials
    Xcat = cat(2, Xtr{:});               % C × (sum T)
    mC   = mean(Xcat,2);
    X0   = double(Xcat) - mC;
    [U,S,~] = svd(X0,'econ');
    s     = diag(S); expl = 100*(s.^2)/sum(s.^2);
    K     = find(cumsum(expl)>=varKeep,1,'first'); if isempty(K), K=size(U,2); end
    Wk    = U(:,1:K);

    Xtr_tem = cellfun(@(x) single(Wk'*(double(x)-mC)), Xtr, 'uni', false);  % K×T
    Xva_tem = cellfun(@(x) single(Wk'*(double(x)-mC)), Xva, 'uni', false);
    inputSize = K;
else
    Xtr_tem = cellfun(@single, Xtr, 'uni', false);    % C×T
    Xva_tem = cellfun(@single, Xva, 'uni', false);
    inputSize = C;
end

Tseq = T;   % fixed length from timelock
end
