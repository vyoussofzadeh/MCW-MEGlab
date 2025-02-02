function [concordanceBA, outlierIndices] = ecpfunc_blandAltmanConcordance(megLI, fmriLI)
% BLANDALTMANCONCORDANCE
%   Returns the proportion of subjects that fall within ±1.96 SD 
%   in a BlandAltman sense (i.e., "concordant" subjects),
%   and optionally which subjects are outliers.

    megLI  = megLI(:);
    fmriLI = fmriLI(:);
    diffVals = megLI - fmriLI;
    nSubj    = length(diffVals);

    meanDiff = mean(diffVals);
    sdDiff   = std(diffVals);
    loaUpper = meanDiff + 1.96 * sdDiff;
    loaLower = meanDiff - 1.96 * sdDiff;

    isOutlier     = (diffVals > loaUpper) | (diffVals < loaLower);
    outlierIndices = find(isOutlier);  % which subjects are outliers

    nOutliers    = sum(isOutlier);
    nInliers     = nSubj - nOutliers;

    % Fraction (or percentage) that are inliers
    concordanceBA = nInliers / nSubj;
end
