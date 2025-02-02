function IQcat = defineIQbins(iqVals)
% Splits IQ into <80 (Below), 80-120 (Average), >120 (Above)
% Adjust thresholds as needed
    IQcatStr = repmat("Average", size(iqVals));
    IQcatStr(iqVals < 80)  = "Below";
    IQcatStr(iqVals > 120) = "Above";
    IQcat = categorical(IQcatStr, ["Below","Average","Above"]);
end
