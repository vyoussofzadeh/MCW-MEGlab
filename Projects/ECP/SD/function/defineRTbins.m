function catArr = defineRTbins(rtArray)
% DEFINERTBINS  Example threshold-based approach for Reaction Times.
% E.g. "Fast" if <1.0s, "Moderate" if 1.0-1.5s, "Slow" if >1.5s.

    catArr = repmat("Slow", size(rtArray));
    catArr(rtArray < 1.0) = "Fast";
    catArr(rtArray >=1.0 & rtArray <=1.5) = "Moderate";
    catArr = categorical(catArr, ["Fast","Moderate","Slow"]);
end
