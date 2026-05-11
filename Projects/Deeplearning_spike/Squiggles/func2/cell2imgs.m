function [X4d, Y1, info] = cell2imgs(Xcell, Y, fixedLen)
%CELL2IMGS Convert {H x T_i} cell epochs to 4-D [H W 1 N] for image models.
%   [X4d, Y1]             = cell2imgs(Xcell, Y)           % auto pad to max T
%   [X4d, Y1, info]       = cell2imgs(Xcell, Y)           % returns sizes
%   [X4d, Y1]             = cell2imgs(Xcell, Y, fixedLen) % pad/crop to fixedLen
%
% Inputs:
%   Xcell   : 1xN cell, each Xcell{i} is [H x T_i] double/single
%   Y       : N x 1 categorical (or convertible)
%   fixedLen: (optional) positive integer. If provided, crops/pads to this.
%
% Outputs:
%   X4d     : [H x W x 1 x N] single
%   Y1      : N x 1 categorical (same order)
%   info    : struct with fields H, W, N, lengths (vector of T_i)

    assert(iscell(Xcell) && ~isempty(Xcell), 'Xcell must be a non-empty cell array.');
    N = numel(Xcell);
    if ~isa(Y,'categorical'); Y = categorical(Y); end
    assert(numel(Y)==N, 'Y must have one label per cell element.');

    % Determine H and T lengths
    H = size(Xcell{1},1);
    T = zeros(N,1);
    for i=1:N
        Xi = Xcell{i};
        assert(size(Xi,1)==H, 'All epochs must have the same #channels (rows).');
        T(i) = size(Xi,2);
    end

    if nargin>=3 && ~isempty(fixedLen)
        W = fixedLen;
    else
        W = max(T);  % auto-pad to longest
    end

    % Allocate 4-D array
    X4d = zeros(H, W, 1, N, 'single');

    % Right-pad (zeros) or crop to W
    for i=1:N
        Xi = single(Xcell{i});
        Ti = size(Xi,2);
        if Ti >= W
            X4d(:,:,1,i) = Xi(:,1:W);
        else
            X4d(:,1:Ti,1,i) = Xi;
        end
    end

    Y1 = Y;

    if nargout>2
        info = struct('H',H,'W',W,'N',N,'lengths',T);
    end
end
