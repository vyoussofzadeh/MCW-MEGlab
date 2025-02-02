function AEDcat = defineAEDbins(aedVals)

AEDcatStr = repmat("3plus", size(aedVals));
AEDcatStr(aedVals == 0) = "0";
AEDcatStr(aedVals == 1) = "1";
AEDcatStr(aedVals == 2) = "2";
AEDcat = categorical(AEDcatStr, ["0","1","2","3plus"]);
end
