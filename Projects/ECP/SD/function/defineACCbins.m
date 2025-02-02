function catArr = defineACCbins(accArray)
% DEFINEACCBINS  Example threshold-based approach for Accuracy.
% E.g. "Low" if <50%, "Mid" if 50-80%, "High" if >80%.

    catArr = repmat("High", size(accArray));
    catArr(accArray < 50) = "Low";
    catArr(accArray >=50 & accArray <=80) = "Mid";
    catArr = categorical(catArr, ["Low","Mid","High"]);
end
