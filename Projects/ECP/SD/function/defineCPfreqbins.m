function CPcat = defineCPfreqbins(cpVals)
% DEFINECPFREQBINS  Convert numeric cp_freq into 4 bins:
%   '0', '1to5', '6to10', '11plus'
%
% EXAMPLE:
%   T1_epil_measures.cp_freq_cat = defineCPfreqbins(T1_epil_measures.CP_freq);

    if ~isvector(cpVals) || ~isnumeric(cpVals)
        error('cpVals must be a numeric vector.');
    end

    CPcatStr = strings(size(cpVals));  % string array for intermediate
    CPcatStr(:) = "0";                 % default
    CPcatStr(cpVals >= 1 & cpVals <= 5)   = "1to5";
    CPcatStr(cpVals >= 6 & cpVals <= 10) = "6to10";
    CPcatStr(cpVals > 10)               = "11plus";

    % Convert to categorical with a defined order
    CPcat = categorical(CPcatStr, ["0","1to5","6to10","11plus"]);
end
