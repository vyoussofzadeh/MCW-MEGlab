function    [trl,trlInfoColDescr, trialinfo, trialSummary, scanStartSamp, scanEndSamp, warninfo] = trialfun_StoryM_neuromag( cfg )
%% This is a wrapper function for the Story/Math core trial definition function trialfun_StoryM_BaseExtractAll.m
% See the help notes of this core function for information on input and
% output variables.

% Copyright (C) 2011-2013 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
%
% This file is part of megconnectome.
%
% megconnectome is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% megconnectome is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with megconnectome.  If not, see <http://www.gnu.org/licenses/>.

% [trl,trialinfo, trlInfoColDescr, trialSummary, scanStartSamp, scanEndSamp, warninfo] = trialfun_StoryM_BaseExtractAll_neuromag( cfg );
[trl,trialinfo, trlInfoColDescr, trialSummary, scanStartSamp, scanEndSamp, warninfo] = trialfun_StoryM_BaseExtractAll_neuromag_corrected( cfg );

end
