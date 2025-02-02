function [kappa, po, pe] = ecpfunc_cohenKappa(confMat)
% cohenKappa  Compute Cohen's kappa for confusion matrix
%   confMat: square confusion matrix (rows = predicted, cols = actual)
%   kappa: computed Cohen's kappa
%   po: observed agreement
%   pe: expected agreement by chance

    total = sum(confMat(:));
    
    % Observed agreement (sum of diagonal / total)
    po = trace(confMat) / total;
    
    % Expected agreement
    rowSums = sum(confMat, 2);  % sum across columns -> row vector
    colSums = sum(confMat, 1);  % sum across rows -> column vector
    pe = (rowSums' * colSums') / (total^2);
    
    % Cohen's kappa
    kappa = (po - pe) / (1 - pe);
end