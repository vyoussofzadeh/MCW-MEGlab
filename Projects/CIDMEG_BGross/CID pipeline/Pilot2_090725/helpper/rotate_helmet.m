function [ChannelMat, R, O] = rotate_helmet(ChannelMat, axisOrR, deg, center, typesToRotate)
% Rotate MEG sensor coil positions (Loc) and orientations (Orient) by a fixed transform.
% Leaves HeadPoints and MRI as-is (i.e., you are rotating the helmet only).
%
% Usage:
%   % 1) In-memory struct:
%   [ChannelMat,R,O] = rotate_helmet(ChannelMat,'z',-90,'lpa_rpa');  % -90° about Z, center at LPA/RPA midpoint
%
%   % 2) On-disk Brainstorm Channel file (fast, no GUI):
%   ChannelFile = '.../data/@raw_myrec/channel_vectorview306.mat';
%   ChannelMat  = in_bst_channel(ChannelFile);
%   [ChannelMat,R,O] = rotate_helmet(ChannelMat,'z',-90,'lpa_rpa');
%   bst_save(ChannelFile, ChannelMat, 'v6');   % save back
%
% Inputs:
%   ChannelMat     : Brainstorm ChannelMat struct (must contain .Channel)
%   axisOrR        : 'x'|'y'|'z'  OR  3x3 rotation matrix
%   deg            : angle in degrees if axisOrR is char (default = -90)
%   center         : 'lpa_rpa' (default) | 'nas' | 'origin' | 'centroid' | [3x1] numeric point (meters)
%   typesToRotate  : cellstr of channel Type(s) to rotate (default = {'MEG','MEG MAG','MEG GRAD'})
%
% Outputs:
%   ChannelMat     : updated struct
%   R              : 3x3 rotation matrix applied
%   O              : 3x1 rotation center used
%
% Notes:
%   - Works for Neuromag/Elekta/MEGIN (multi-coil sensors: Loc is 3xN per channel).
%   - Recompute the head model after saving the updated Channel file.

if nargin < 2 || isempty(axisOrR), axisOrR = 'z'; end
if nargin < 3 || isempty(deg),     deg     = -90; end
if nargin < 4 || isempty(center),  center  = 'lpa_rpa'; end
if nargin < 5 || isempty(typesToRotate), typesToRotate = {'MEG','MEG MAG','MEG GRAD'}; end

assert(isfield(ChannelMat,'Channel') && ~isempty(ChannelMat.Channel), 'ChannelMat.Channel missing');

% -------- Build rotation matrix R --------
if isnumeric(axisOrR) && isequal(size(axisOrR),[3,3])
    R = double(axisOrR);
else
    th = deg*pi/180;
    switch lower(axisOrR)
        case 'x', R = [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];
        case 'y', R = [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
        case 'z', R = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
        otherwise, error('axisOrR must be ''x'',''y'',''z'' or a 3x3 matrix');
    end
end

% -------- Choose rotation center O --------
O = [];
% Try LPA/RPA/NAS from HeadPoints if available
if isfield(ChannelMat,'HeadPoints') && ~isempty(ChannelMat.HeadPoints) ...
        && isfield(ChannelMat.HeadPoints,'Loc')
    hp = ChannelMat.HeadPoints;
    Lhp = double(hp.Loc);
    switch lower(center)
        case 'lpa_rpa'
            iL = find(strcmp(hp.Label,'LPA'),1);
            iR = find(strcmp(hp.Label,'RPA'),1);
            if ~isempty(iL) && ~isempty(iR), O = 0.5*(Lhp(:,iL) + Lhp(:,iR)); end
        case 'nas'
            iN = find(strcmp(hp.Label,'NAS'),1);
            if ~isempty(iN), O = Lhp(:,iN); end
    end
end
if isempty(O)
    switch lower(center)
        case 'origin',   O = [0;0;0];
        case 'centroid', O = local_centroid(ChannelMat, typesToRotate);
        otherwise
            if isnumeric(center) && isequal(size(center),[3,1])
                O = double(center);
            else
                % default fallback if LPA/RPA not found:
                O = local_centroid(ChannelMat, typesToRotate);
            end
    end
end

% -------- Rotate selected channel types --------
isTargetType = @(tp) any(strcmpi(strtrim(tp), typesToRotate));
for k = 1:numel(ChannelMat.Channel)
    ch = ChannelMat.Channel(k);
    if isempty(ch.Loc) || ~isTargetType(ch.Type), continue; end

    L = double(ch.Loc);                         % 3 x nCoils
    L = R*(L - O) + O;                          % rotate about O
    ChannelMat.Channel(k).Loc = cast(L, 'like', ch.Loc);

    if ~isempty(ch.Orient)
        Q = double(ch.Orient);
        ChannelMat.Channel(k).Orient = cast(R*Q, 'like', ch.Orient);
    end
end

% -------- Tag the comment (optional) --------
try
    if isfield(ChannelMat,'Comment') && ischar(ChannelMat.Comment)
        ChannelMat.Comment = sprintf('%s | rot-helmet (%s,%g°)', ChannelMat.Comment, ...
                                     ischar(axisOrR)*axisOrR + ~ischar(axisOrR)*'R', deg);
    end
end

end % rotate_helmet

% ---------- helpers ----------
function C = local_centroid(ChannelMat, typesToRotate)
    C = [0;0;0]; n = 0;
    isTargetType = @(tp) any(strcmpi(strtrim(tp), typesToRotate));
    for kk = 1:numel(ChannelMat.Channel)
        ch = ChannelMat.Channel(kk);
        if isempty(ch.Loc) || ~isTargetType(ch.Type), continue; end
        L = double(ch.Loc);
        C = C + sum(L,2);
        n = n + size(L,2);
    end
    if n>0, C = C / n; else, C = [0;0;0]; end
end
