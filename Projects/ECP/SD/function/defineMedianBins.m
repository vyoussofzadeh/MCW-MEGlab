function catArr = defineMedianBins(numericArr)
% DEFINEMEDIANBINS  Splits a numeric array into two categories:
%   "<=median" or ">median"
% Returns a categorical array of the same size.

    medVal = nanmedian(numericArr);
    catArr = repmat("<=Median", size(numericArr));
    catArr(numericArr > medVal) = ">Median";
    catArr = categorical(catArr, ["<=Median", ">Median"]);
end
