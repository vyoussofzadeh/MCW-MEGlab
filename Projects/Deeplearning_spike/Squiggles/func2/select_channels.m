function idxKeep = select_channels(scores, K, chLabels)
% scores: Cx1 vector; K: number to keep; chLabels (optional) cellstr
C = numel(scores);
K = min(K, C);
if nargin<3 || isempty(chLabels)
    [~,ord] = sort(scores,'descend');
    idxKeep = sort(ord(1:K)); return;
end

isMag  = endsWith(chLabels,'1');
isGrad = endsWith(chLabels,'2') | endsWith(chLabels,'3');

% Proportional split
nMag  = sum(isMag);  nGrad = sum(isGrad);
kMag  = round(K * nMag / max(nMag+nGrad,1));
kGrad = K - kMag;

idxMag  = find(isMag);
idxGrad = find(isGrad);

[~,ordM] = sort(scores(idxMag),'descend');
[~,ordG] = sort(scores(idxGrad),'descend');

keepM = idxMag(ordM(1:min(kMag, numel(ordM))));
keepG = idxGrad(ordG(1:min(kGrad, numel(ordG))));

idxKeep = sort([keepM; keepG]);
% if we came up short (e.g., missing labels), fill from global top
if numel(idxKeep) < K
    [~,ord] = sort(scores,'descend');
    fill = setdiff(ord, idxKeep, 'stable');
    idxKeep = sort([idxKeep; fill(1:(K-numel(idxKeep)))]);
end
end
