function EHQcat = defineHandedness(ehqVals, threshold)
% DEFINEHANDEDNESS  Converts numeric EHQ scores (-100..+100)
%   into categorical Left/Ambi/Right based on a threshold.
%
% USAGE:
%   EHQcat = defineHandedness(ehqVals, 40);
%
% INPUT:
%   ehqVals  : Nx1 numeric array (Edinburgh Handedness)
%   threshold: scalar (e.g., 40). 
%              If EHQ > +threshold => "Right"
%              If EHQ < -threshold => "Left"
%              Else => "Ambi"
%
% OUTPUT:
%   EHQcat   : Nx1 categorical array, categories {'Left','Ambi','Right'}
%
% Example:
%   T1_epil_measures.EHQcat = defineHandedness(T1_epil_measures.EHQ, 40);

    if ~isvector(ehqVals) || ~isnumeric(ehqVals)
        error('ehqVals must be a numeric vector.');
    end

    % Start with "Ambi" for all
    EHQcatStr = repmat("Ambi", size(ehqVals));

    % Assign "Right" where EHQ > threshold
    EHQcatStr(ehqVals > threshold) = "Right";

    % Assign "Left" where EHQ < -threshold
    EHQcatStr(ehqVals < -threshold) = "Left";

    % Convert to categorical with a set ordering
    EHQcat = categorical(EHQcatStr, ["Left","Ambi","Right"]);

end
