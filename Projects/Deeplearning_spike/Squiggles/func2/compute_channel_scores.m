function [scores, perChanFeat] = compute_channel_scores(Xcells, Y, method)
% Compute an informativeness score per channel using TRAIN data.
% Xcells: 1xN cell, each [C x T] (normalized epochs)
% Y     : N x 1 categorical (e.g., {'NoSpike','Spike'})
% method: 'fisher' (default), 'rms', 'mi'

if nargin<3 || isempty(method), method = 'fisher'; end
assert(iscell(Xcells)&&~isempty(Xcells),'Xcells must be a nonempty cell.');
C = size(Xcells{1},1);
N = numel(Xcells);

% Reduce each epoch to a per-channel scalar (fast, robust):
% use RMS across time per channel (other choices: mean|x|, peak2peak)
perChanFeat = zeros(N, C, 'single');
for i = 1:N
    Xi = single(Xcells{i});                  % [C x T]
    perChanFeat(i,:) = rms(Xi, 2);           % 1 x C
end

switch lower(method)
    case 'fisher'
        % Fisher score per channel: (mu1 - mu0)^2 / (var1 + var0 + eps)
        cats = categories(Y);
        assert(numel(cats)==2, 'Fisher scoring assumes 2 classes.');
        m0 = mean(perChanFeat(Y==cats{1},:), 1);
        m1 = mean(perChanFeat(Y==cats{2},:), 1);
        v0 = var( perChanFeat(Y==cats{1},:), 0, 1);
        v1 = var( perChanFeat(Y==cats{2},:), 0, 1);
        scores = ((m1 - m0).^2) ./ (v1 + v0 + 1e-12);
    case 'rms'
        % Larger RMS across all epochs = "louder" channels (useful if pre-denoised)
        scores = mean(perChanFeat, 1);
    case 'mi'
        % Mutual information between per-channel feature and label (coarse bins)
        cats = categories(Y); ybin = double(Y==cats{2});  % 0/1
        scores = zeros(1,C);
        for c = 1:C
            x = perChanFeat(:,c);
            edges = prctile(x,[0 20 40 60 80 100]); edges(1) = edges(1)-eps;
            xb = discretize(x, unique(edges));
            scores(c) = mutual_info_discrete(xb, ybin);
        end
    otherwise
        error('Unknown method: %s', method);
end

scores = double(scores(:));  % C x 1
end

function I = mutual_info_discrete(xb, y)
% xb: discrete 1..B, y: 0/1
B = max(xb); I = 0;
px = accumarray(xb,1,[B 1])/numel(xb);
py = [mean(y==0); mean(y==1)];
for b = 1:B
    idx = (xb==b);
    pxy = [mean(y(idx)==0); mean(y(idx)==1)] * px(b);
    pxpy = px(b) * py;
    mask = pxy>0;
    I = I + sum(pxy(mask) .* log(pxy(mask)./pxpy(mask)));
end
I = max(I,0);
end
