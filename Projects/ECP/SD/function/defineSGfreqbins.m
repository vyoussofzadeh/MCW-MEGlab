function SGcat = defineSGfreqbins(sgVals)
    if ~isvector(sgVals) || ~isnumeric(sgVals)
        error('sgVals must be numeric.');
    end

    SGcatStr = repmat("0", size(sgVals));
    SGcatStr(sgVals >= 1 & sgVals <= 2) = "1to2";
    SGcatStr(sgVals > 2) = "3plus";

    SGcat = categorical(SGcatStr, ["0","1to2","3plus"]);
end
