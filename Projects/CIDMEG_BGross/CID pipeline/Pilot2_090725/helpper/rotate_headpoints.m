% Rotate Brainstorm head points by a fixed transform, and (optionally) rotate
% Channel coil positions/orientations by the same transform to keep geometry consistent.
%
% Usage:
%   neuromag_channs = rotate_headpoints(neuromag_channs)                          % default: -90 deg about Z, EXTRA only
%   neuromag_channs = rotate_headpoints(neuromag_channs,'z',-90,'ALL',true)       % rotate ALL head points + Channels
%   [neuromag_channs,R,O] = rotate_headpoints(neuromag_channs,eye(3),[],[],true)  % pass a custom 3x3 R; rotate Channels too
%
% Inputs:
%   neuromag_channs  : Brainstorm ChannelMat-like struct with HeadPoints (and optionally Channel)
%   axisOrR          : 'x'|'y'|'z'  OR a 3x3 numeric rotation matrix (default 'z')
%   deg              : rotation in degrees if axisOrR is a char axis (default -90). Ignored if axisOrR is 3x3.
%   rotateWhich      : 'EXTRA' (default) | 'ALL' | logical(1xN) mask over head points
%   rotateChannelsToo: logical, also apply to Channel(k).Loc / .Orient (default false)
%
% Outputs:
%   neuromag_channs  : struct with rotated HeadPoints (and Channels if requested)
%   R                : 3x3 rotation matrix applied
%   O                : 3x1 rotation center used (midpoint of LPA/RPA when available; else centroid of selected points)
%
function [neuromag_channs, R, O] = rotate_headpoints(neuromag_channs, axisOrR, deg, rotateWhich, rotateChannelsToo)

% ---- defaults ----
if nargin < 2 || isempty(axisOrR),        axisOrR = 'z';   end
if nargin < 3 || isempty(deg),            deg      = -90;  end
if nargin < 4 || isempty(rotateWhich),    rotateWhich = 'EXTRA'; end
if nargin < 5 || isempty(rotateChannelsToo), rotateChannelsToo = false; end

% ---- sanity checks ----
if ~isfield(neuromag_channs,'HeadPoints') || isempty(neuromag_channs.HeadPoints) ...
        || ~isfield(neuromag_channs.HeadPoints,'Loc')
    error('rotate_headpoints:MissingHeadPoints','No HeadPoints found in input struct.');
end

hp = neuromag_channs.HeadPoints;
L  = double(hp.Loc);                       % 3 x N
N  = size(L,2);
if N ~= numel(hp.Type)
    error('rotate_headpoints:SizeMismatch','HeadPoints.Loc and HeadPoints.Type size mismatch.');
end

% ---- pick indices to rotate ----
if ischar(rotateWhich) || isstring(rotateWhich)
    switch upper(string(rotateWhich))
        case "EXTRA", idx = strcmpi(hp.Type,'EXTRA');
        case "ALL",   idx = true(1,N);
        otherwise, error('rotate_headpoints:BadRotateWhich','rotateWhich must be ''EXTRA'', ''ALL'', or a logical mask.');
    end
elseif islogical(rotateWhich) && isvector(rotateWhich) && numel(rotateWhich)==N
    idx = reshape(rotateWhich,1,[]);
else
    error('rotate_headpoints:BadRotateWhich','rotateWhich must be ''EXTRA'', ''ALL'', or a logical mask of length N.');
end

% ---- rotation center O: midpoint(LPA,RPA) if present; else centroid of selected points ----
iL = find(strcmp(hp.Label,'LPA'), 1);
iR = find(strcmp(hp.Label,'RPA'), 1);
if ~isempty(iL) && ~isempty(iR)
    O = 0.5*(L(:,iL) + L(:,iR));          % 3x1 origin
else
    O = mean(L(:,idx),2);                  % centroid of selected points
end

% ---- rotation matrix R ----
if isnumeric(axisOrR) && isequal(size(axisOrR),[3 3])
    R = double(axisOrR);
else
    th = deg*pi/180;                       % degrees->radians
    switch lower(char(axisOrR))
        case 'x', R = [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];
        case 'y', R = [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
        case 'z', R = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
        otherwise, error('rotate_headpoints:BadAxis','axisOrR must be ''x'',''y'',''z'' or a 3x3 numeric matrix.');
    end
end

% ---- apply rotation to selected head points about O ----
L(:,idx) = R*(L(:,idx) - O) + O;

% ---- write back (preserve numeric class) ----
hp.Loc = cast(L, 'like', neuromag_channs.HeadPoints.Loc);
neuromag_channs.HeadPoints = hp;

% ---- optional: rotate Channel coil positions and orientations by same R/O ----
if rotateChannelsToo && isfield(neuromag_channs,'Channel') && ~isempty(neuromag_channs.Channel)
    for k = 1:numel(neuromag_channs.Channel)
        % Positions (each column is a coil location)
        Lk = neuromag_channs.Channel(k).Loc;
        if ~isempty(Lk)
            Lkd = double(Lk);
            LkR = R*(Lkd - O) + O;
            neuromag_channs.Channel(k).Loc = cast(LkR, 'like', Lk);
        end
        % Orientations (no translation)
        Ok = neuromag_channs.Channel(k).Orient;
        if ~isempty(Ok)
            Okd = double(Ok);
            neuromag_channs.Channel(k).Orient = cast(R*Okd, 'like', Ok);
        end
    end
end
end


% % Rotate Brainstorm head points by a fixed angle.
% % - Rotates only EXTRA (scalp) points by default.
% % - Rotation center: midpoint of LPA/RPA (head origin); falls back to
% %   centroid of the chosen points if LPA/RPA are missing.
% 
% function neuromag_channs = rotate_headpoints(neuromag_channs, axisChar, deg, rotateWhich)
% 
% if nargin < 4, rotateWhich = 'EXTRA'; end  % 'EXTRA' | 'ALL'
% if nargin < 3, deg = -90; end              % default -90 degrees
% if nargin < 2, axisChar = 'z'; end         % 'x' | 'y' | 'z'
% 
% hp = neuromag_channs.HeadPoints;
% L  = double(hp.Loc);                        % 3 x N
% N  = numel(hp.Type);
% 
% % --- pick indices to rotate ---
% switch upper(rotateWhich)
%   case 'EXTRA', idx = strcmpi(hp.Type,'EXTRA');
%   case 'ALL',   idx = true(1,N);
%   otherwise, error('rotateWhich must be ''EXTRA'' or ''ALL''.');
% end
% 
% % --- rotation center: midpoint(LPA,RPA) if present; else centroid ---
% iL = find(strcmp(hp.Label,'LPA'), 1);
% iR = find(strcmp(hp.Label,'RPA'), 1);
% if ~isempty(iL) && ~isempty(iR)
%     O = 0.5*(L(:,iL) + L(:,iR));           % 3x1 origin
% else
%     O = mean(L(:,idx),2);                   % centroid of selected points
% end
% 
% % --- rotation matrix (right-handed, degrees->radians) ---
% th = deg*pi/180;
% switch lower(axisChar)
%   case 'x'
%     R = [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];
%   case 'y'
%     R = [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
%   case 'z'
%     R = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
%   otherwise
%     error('axisChar must be ''x'',''y'', or ''z''.');
% end
% 
% % --- apply rotation about O to selected points ---
% L(:,idx) = R*(L(:,idx) - O) + O;
% 
% % --- write back (preserve single precision) ---
% hp.Loc = single(L);
% neuromag_channs.HeadPoints = hp;
% 
% % ---- OPTIONAL: rotate sensor locations/orientations by the same R ----
% % If this struct is a full ChannelMat and you want sensors to follow:
% % if isfield(neuromag_channs,'Channel')
% %   for k = 1:numel(neuromag_channs.Channel)
% %     if ~isempty(neuromag_channs.Channel(k).Loc)
% %       neuromag_channs.Channel(k).Loc = single(R*(double(neuromag_channs.Channel(k).Loc) - O) + O);
% %     end
% %     if ~isempty(neuromag_channs.Channel(k).Orient)
% %       neuromag_channs.Channel(k).Orient = single(R*double(neuromag_channs.Channel(k).Orient));
% %     end
% %   end
% % end
% end
