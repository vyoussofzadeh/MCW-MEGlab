function Xout = maggrad_median_scale_per_epoch(Xcells, refLbl)
% Scale each epoch so mags/grads are on comparable ranges using per-epoch medians
isMag  = endsWith(refLbl,'1');
isGrad = endsWith(refLbl,'2') | endsWith(refLbl,'3');

Xout = cell(size(Xcells));
for i = 1:numel(Xcells)
    X = Xcells{i};
    mMag  = median(abs(X(isMag ,:)), 'all','omitnan') + eps;
    mGrad = median(abs(X(isGrad,:)), 'all','omitnan') + eps;
    X(isMag ,:) = X(isMag ,:) ./ mMag;
    X(isGrad,:) = X(isGrad,:) ./ mGrad;
    Xout{i} = X;
end
end
