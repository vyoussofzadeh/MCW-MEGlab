function dlnet = pretrain_masked_recon_resnet(trainImgs, varargin)
% Self-supervised masked reconstruction pretrain (MATLAB-version friendly).
% trainImgs: [C x W x 1 x N] single
% Params: 'MaskFrac',0.2, 'Epochs',20, 'MB',128, 'LR',3e-4, 'Beta1',0.9, 'Beta2',0.999, 'Epsilon',1e-8

p = inputParser;
addParameter(p,'MaskFrac',0.2);
addParameter(p,'Epochs',20);
addParameter(p,'MB',128);
addParameter(p,'LR',3e-4);
addParameter(p,'Beta1',0.9);
addParameter(p,'Beta2',0.999);
addParameter(p,'Epsilon',1e-8);
parse(p,varargin{:});
mf = p.Results.MaskFrac; E = p.Results.Epochs; MB = p.Results.MB; 
LR = p.Results.LR; beta1 = p.Results.Beta1; beta2 = p.Results.Beta2; eps0 = p.Results.Epsilon;

[C,W,~,N] = size(trainImgs);

% Backbone (ensure you've updated emsnet_lite_backbone.m as we discussed)
lgraph = emsnet_lite_backbone(C,W);
dlnet  = dlnetwork(lgraph);

% Initialize Adam state (same shape as Learnables)
[mState, vState, t] = initAdamState(dlnet);

numIter = ceil(N/MB);

for epoch = 1:E
    idx = randperm(N);
    for it = 1:numIter
        batch = idx((it-1)*MB+1 : min(it*MB,N));
        X = trainImgs(:,:,:,batch);                     % [C W 1 B]
        [Xmask, target, M] = make_mask_batch(X, mf);    % mask columns

        dlX = dlarray(Xmask,'SSCB');
        dlY = dlarray(target,'SSCB');
        dlM = dlarray(M,'SSCB');

        [loss, grads] = dlfeval(@maskedReconLoss, dlnet, dlX, dlY, dlM);

        % Custom Adam step (version-agnostic)
        t = t + 1;
        [dlnet, mState, vState] = adam_step(dlnet, grads, mState, vState, t, LR, beta1, beta2, eps0);
    end
    fprintf('Pretrain %d/%d  loss=%.5f\n', epoch, E, double(gather(extractdata(loss))));
end
end

% ==== helpers ====

function [loss, grads] = maskedReconLoss(dlnet, dlX, dlY, dlM)
Yhat = forward(dlnet, dlX);           % expects output layer named 'feat_out'
E = (Yhat - dlY) .* dlM;              % only masked positions
numMask = sum(dlM(:)) + eps;
loss = sum(E(:).^2) / numMask;
grads = dlgradient(loss, dlnet.Learnables);
end

function [Xmask, target, M] = make_mask_batch(X, maskFrac)
% X: [C W 1 B]; maskFrac in (0,1). Masks whole time columns (broadcast across channels).
[C,W,~,B] = size(X);
target = X;
M      = zeros(1,W,1,B,'like',X);
Xmask  = X;
for b=1:B
    k = max(1, round(maskFrac * W));
    cols = randperm(W, k);
    M(1,cols,1,b) = 1;
    Xmask(:,cols,1,b) = 0;
end
end

function [mState, vState, t] = initAdamState(dlnet)
L = dlnet.Learnables.Value;
mState = cell(size(L));
vState = cell(size(L));
for i=1:numel(L)
    mState{i} = zeros(size(L{i}), 'like', L{i});
    vState{i} = zeros(size(L{i}), 'like', L{i});
end
t = 0;
end

function [dlnet, mState, vState] = adam_step(dlnet, grads, mState, vState, t, lr, b1, b2, eps0)
% Version-agnostic Adam update for dlnetwork
for i = 1:numel(dlnet.Learnables.Value)
    g = grads.Value{i};
    if isempty(g), continue; end                % some layers may have no grads
    m = mState{i}; v = vState{i};

    m = b1*m + (1-b1)*g;
    v = b2*v + (1-b2)*(g.^2);

    mhat = m ./ (1 - b1^t);
    vhat = v ./ (1 - b2^t);

    param = dlnet.Learnables.Value{i};
    param = param - lr * mhat ./ (sqrt(vhat) + eps0);

    dlnet.Learnables.Value{i} = param;
    mState{i} = m; vState{i} = v;
end
end
