function Xw = warp_sequences_to_length_interp(Xcells, Tfix, method)
% Interpolate each C×T epoch onto a uniform Tfix grid in [0,1].
% method: 'linear' (default) | 'pchip' | 'spline' | 'nearest'
if nargin < 3, method = 'linear'; end
tnew = linspace(0,1,Tfix);
Xw = cell(size(Xcells));
for i = 1:numel(Xcells)
    x = Xcells{i};                 % C×T
    Ti = size(x,2);
    if Ti == Tfix
        Xw{i} = x;
        continue;
    end
    told = linspace(0,1,Ti);
    % vectorized over channels via transpose
    Xw{i} = interp1(told, double(x).', tnew, method, 'extrap').';  % C×Tfix
end
end



% function Xw = warp_sequences_to_length_interp(Xcells, Tfix, method)
% if nargin<3, method='pchip'; end
% tnew = linspace(0,1,Tfix);
% Xw = cell(size(Xcells));
% for i=1:numel(Xcells)
%   x = Xcells{i}; Ti = size(x,2);
%   if Ti==Tfix, Xw{i}=x; else
%     told = linspace(0,1,Ti);
%     Xw{i} = interp1(told, double(x).', tnew, method, 'extrap').';
%   end
% end
% end
