function net = run_A3_pretrain_and_finetune(trainImgs, trainLabs, valImgs, valLabs, varargin)
% RUN_A3_PRETRAIN_AND_FINETUNE
% Self-supervised masked-reconstruction pretraining (A3) + supervised fine-tuning.
%
% Inputs:
%   trainImgs : [C x W x 1 x Ntr] single
%   trainLabs : Ntr x 1 categorical
%   valImgs   : [C x W x 1 x Nv]  single
%   valLabs   : Nv x 1 categorical
%
% Name-Value:
%   'MaskFrac'        (0.2)   fraction of time columns to mask during pretrain
%   'PretrainEpochs'  (10)
%   'FT_Epochs'       (30)
%   'MB'              (128)   minibatch size
%   'LR'              (3e-4)  learning rate
%   'UseClassWeights' (true)
%
% Output:
%   net : trained classification DAGNetwork

% ----------- Parse args -----------
p = inputParser;
addParameter(p,'MaskFrac',0.2);
addParameter(p,'PretrainEpochs',10);
addParameter(p,'FT_Epochs',30);
addParameter(p,'MB',128);
addParameter(p,'LR',3e-4);
addParameter(p,'UseClassWeights',true);
parse(p,varargin{:});
mf   = p.Results.MaskFrac;
Epre = p.Results.PretrainEpochs;
Eft  = p.Results.FT_Epochs;
MB   = p.Results.MB;
LR   = p.Results.LR;
useW = p.Results.UseClassWeights;

% ----------- Sanity -----------
assert(ndims(trainImgs)==4 && size(trainImgs,3)==1, 'Expect trainImgs as [C x W x 1 x N].');
assert(isa(trainImgs,'single') && isa(valImgs,'single'), 'Use single precision arrays.');
cats = categories(trainLabs);
trainLabs = reordercats(trainLabs, cats);
valLabs   = reordercats(valLabs,   cats);

C = size(trainImgs,1);  W = size(trainImgs,2);

% ----------- A3.1: Build backbone -----------
lgraph = emsnet_lite_backbone(C,W);
dlnet  = dlnetwork(lgraph);

% ----------- A3.2: Self-supervised pretraining (masked recon) -----------
dlnet = pretrain_masked_recon(dlnet, trainImgs, mf, Epre, MB, LR);

% ----------- A3.3: Attach classifier head & (optionally) class weights -----------
lgraph = layerGraph(dlnet);  % convert current net to layerGraph
cls = [
    globalAveragePooling2dLayer("Name","gap")
    dropoutLayer(0.3,"Name","drop")
    fullyConnectedLayer(numel(cats),"Name","fc")
    softmaxLayer("Name","sm")
    classificationLayer("Name","cls")
];
lgraph = addLayers(lgraph, cls);
lgraph = connectLayers(lgraph, "feat_out", "gap");

if useW
    cnt = countcats(trainLabs);
    w   = sum(cnt)./max(cnt,1); w = w/mean(w);
    clsW = classificationLayer('Name','cls','Classes',cats,'ClassWeights',w);
    lgraph = replaceLayer(lgraph, "cls", clsW);
end

% ----------- A3.4: Fine-tune -----------
valFreq = max(10, ceil(size(trainImgs,4)/MB));
opts = trainingOptions('adam', ...
    'InitialLearnRate',LR, ...
    'MaxEpochs',Eft, ...
    'MiniBatchSize',MB, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{valImgs, valLabs}, ...
    'ValidationFrequency',valFreq, ...
    'ValidationPatience',10, ...
    'ExecutionEnvironment','auto', ...
    'Plots','training-progress', ...
    'Verbose',true);

net = trainNetwork(trainImgs, trainLabs, lgraph, opts);

end % main


% ========== Helpers below (same file) ==========

function lgraph = emsnet_lite_backbone(C,W)
% Input: imageInputLayer([C W 1]); Output tap: layer "feat_out"
lgraph = layerGraph();
lgraph = addLayers(lgraph, imageInputLayer([C W 1], "Normalization","none", "Name","in"));

% 2D spatio-temporal path
g = [
    convolution2dLayer([7 7], 32, "Padding","same", "Name","g1")
    reluLayer("Name","r1")
    convolution2dLayer([3 3], 64, "Padding","same", "Name","g2")
    reluLayer("Name","r2")
];
lgraph = addLayers(lgraph, g);
lgraph = connectLayers(lgraph, "in", "g1");

% 1D temporal path (shared across sensors)
t = [
    convolution2dLayer([1 7], 32, "Padding","same", "Name","t1")
    reluLayer("Name","tr1")
    convolution2dLayer([1 5], 64, "Padding","same", "Name","t2")
    reluLayer("Name","tr2")
];
lgraph = addLayers(lgraph, t);
lgraph = connectLayers(lgraph, "in", "t1");

% fuse both paths
lgraph = addLayers(lgraph, depthConcatenationLayer(2, "Name","cat"));
lgraph = connectLayers(lgraph, "r2",  "cat/in1");
lgraph = connectLayers(lgraph, "tr2", "cat/in2");

% projection head; keep spatial size; output feature map "feat_out"
h = [
    convolution2dLayer([1 1], 64, "Padding","same", "Name","proj")
    reluLayer("Name","hr")
    convolution2dLayer([3 3], 32, "Padding","same", "Name","up1")
    reluLayer("Name","hr2")
    convolution2dLayer([1 1], 1,  "Padding","same", "Name","feat_out")
];
lgraph = addLayers(lgraph, h);
lgraph = connectLayers(lgraph, "cat", "proj");
end


function dlnet = pretrain_masked_recon(dlnet, trainImgs, maskFrac, E, MB, LR)
% Masked reconstruction pretraining on trainImgs.
[C,W,~,N] = size(trainImgs);
numIter = ceil(N/MB);
vel=[]; ag=[]; asq=[];

for epoch = 1:E
    idx = randperm(N);
    for it=1:numIter
        b = idx((it-1)*MB+1:min(it*MB,N));
        X = trainImgs(:,:,:,b);                        % [C W 1 B]
        [Xmask, target, M] = make_mask_batch(X, maskFrac);

        dlX = dlarray(Xmask,'SSCB');
        dlY = dlarray(target,'SSCB');
        dlM = dlarray(M,'SSCB');

        [loss, grads] = dlfeval(@maskedReconLoss, dlnet, dlX, dlY, dlM);
        [dlnet,vel,ag,asq] = adamupdate(dlnet, grads, vel, ag, asq, ...
                                        (epoch-1)*numIter+it, LR, 0.9, 0.999, 1e-8);
    end
    fprintf('Pretrain %d/%d  loss=%.4f\n', epoch, E, double(gather(extractdata(loss))));
end
end


function [loss, grads] = maskedReconLoss(dlnet, dlX, dlY, dlM)
Yhat = forward(dlnet, dlX);            % [C W 1 B], from "feat_out"
E = (Yhat - dlY) .* dlM;               % only masked columns contribute
numMask = sum(dlM(:)) + eps;
loss = sum(E(:).^2) / numMask;         % masked MSE
grads = dlgradient(loss, dlnet.Learnables);
end


function [Xmask, target, M] = make_mask_batch(X, maskFrac)
% X: [C W 1 B]; masks whole time-columns per sample
[C,W,~,B] = size(X);
target = X;
M = zeros(1,W,1,B,'like',X);           % broadcast across C
Xmask = X;
k = max(1, round(maskFrac * W));
for b=1:B
    cols = randperm(W, k);
    M(1,cols,1,b) = 1;
    Xmask(:,cols,1,b) = 0;
end
end
