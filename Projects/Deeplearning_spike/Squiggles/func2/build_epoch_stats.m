function F = build_epoch_stats(Xcells, varargin)
% BUILD_EPOCH_STATS  Turn {C×T} epochs into an N×D feature matrix.
% Default features per channel: mean, std, peak2peak  ? D = 3*C
% Optional: add beta bandpower (set 'fs' to enable)
%
% F = build_epoch_stats(Xcells, 'fs', fs, 'addBeta', true)

p = inputParser;
addParameter(p,'fs',[],@(x) isempty(x) || (isscalar(x)&&x>0));
addParameter(p,'addBeta',false,@islogical);
parse(p,varargin{:});
fs = p.Results.fs; addBeta = p.Results.addBeta;

N = numel(Xcells);
C = size(Xcells{1},1);

perChan = @(X) [ ...
    mean(X,2,'omitnan'), ...
    std(X,0,2,'omitnan'), ...
    (max(X,[],2) - min(X,[],2)) ...
    ]; % -> C×3

if addBeta && ~isempty(fs)
    % simple bandpower 1330 Hz using periodogram
    getBP = @(x) bandpower(x', fs, [13 30])'; % returns C×1
else
    getBP = [];
end

% preallocate (worst-case with beta)
D = 3*C + (addBeta && ~isempty(fs))*C;
F = zeros(N, D, 'single');

for i = 1:N
    Xi = single(Xcells{i});         % C×T
    A = perChan(Xi);                % C×3
    if ~isempty(getBP)
        bp = getBP(Xi);             % C×1
        Fi = [A, bp];               % C×4
    else
        Fi = A;                     % C×3
    end
    F(i,:) = reshape(Fi, 1, []);    % 1×(C*featuresPerChan)
end

% clean NaN/Inf (rare, e.g., empty windows)
F(~isfinite(F)) = 0;
end
