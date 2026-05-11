function [net, dlnet_pre, cm, acc, bal] = run_A3_finetune_masked(XTrain2, YTrain, XVal2, YVal, varargin)
% RUN_A3_FINETUNE_MASKED  Self-supervised pretrain (masked recon) + fine-tune classifier.
% Inputs:
%   XTrain2, XVal2 : cell of [C x T] (normalized epochs)
%   YTrain,  YVal  : categorical labels with {'NoSpike','Spike'} etc.
% Params (name/value):
%   'FixedLenPct'  : percentile of train lengths for packing (default 85)
%   'MaskFrac'     : fraction of time columns masked in pretrain (default 0.2)
%   'PreEpochs'    : pretrain epochs (default 10)
%   'FT_Epochs'    : fine-tune epochs (default 30)
%   'MB'           : minibatch size (default 128)
%   'LR'           : learning rate (default 3e-4)

p = inputParser;
addParameter(p,'FixedLenPct',85);
addParameter(p,'MaskFrac',0.20);
addParameter(p,'PreEpochs',10);
addParameter(p,'FT_Epochs',30);
addParameter(p,'MB',128);
addParameter(p,'LR',3e-4);
parse(p,varargin{:});
pct = p.Results.FixedLenPct; mf = p.Results.MaskFrac;
Epre = p.Results.PreEpochs; Eft = p.Results.FT_Epochs;
MB = p.Results.MB; LR = p.Results.LR;

%% 1) Pack cells -> [C x W x 1 x N]
fixedLen = round(prctile(cellfun(@(x) size(x,2), XTrain2), pct));
[trainImgs, trainLabs] = cell2imgs(XTrain2, YTrain, fixedLen);
[valImgs,   valLabs]   = cell2imgs(XVal2,   YVal,   fixedLen);

% Keep class order consistent
cats = categories(trainLabs);
trainLabs = reordercats(trainLabs, cats);
valLabs   = reordercats(valLabs,   cats);

C = size(trainImgs,1); W = size(trainImgs,2);

%% 2) Self-supervised pretraining (masked reconstruction)
dlnet_pre = pretrain_masked_recon_resnet(trainImgs, ...
    'MaskFrac', mf, 'Epochs', Epre, 'MB', MB, 'LR', LR);

%% 3) Attach classifier head and fine-tune
% Convert to layerGraph, add GAP?Dropout?FC?Softmax?ClassLayer
lgraph = layerGraph(dlnet_pre);

% Class weights from TRAIN labels
cnt = countcats(trainLabs);
w   = sum(cnt)./max(cnt,1); w = w/mean(w);
clsLayer = classificationLayer('Name','cls','Classes',cats,'ClassWeights',w);

cls = [
    globalAveragePooling2dLayer("Name","gap")
    dropoutLayer(0.3,"Name","drop")
    fullyConnectedLayer(numel(cats),"Name","fc")
    softmaxLayer("Name","sm")
    classificationLayer("Name","cls")];
lgraph = addLayers(lgraph, cls);
lgraph = replaceLayer(lgraph, "cls", clsLayer);          % weighted class layer
lgraph = connectLayers(lgraph, "feat_out", "gap");       % tap from backbone

% Training options
valFreq = max(10, ceil(size(trainImgs,4)/MB));
opts = trainingOptions('adam', ...
    'InitialLearnRate',LR, 'MaxEpochs',Eft, 'MiniBatchSize',MB, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{valImgs, valLabs}, ...
    'ValidationFrequency', valFreq, 'ValidationPatience', 10, ...
    'ExecutionEnvironment','auto', ...
    'Plots','training-progress', 'Verbose', true);

% Fine-tune
net = trainNetwork(trainImgs, trainLabs, lgraph, opts);

%% 4) Validation metrics
Yp = classify(net, valImgs);
% unify types/order
valLabsU = categorical(valLabs, cats, 'Ordinal', false);
YpU      = categorical(Yp,      cats, 'Ordinal', false);
cm  = confusionmat(valLabsU, YpU, 'Order', cats);
acc = sum(diag(cm))/sum(cm,'all');
rec = diag(cm) ./ max(1, sum(cm,2));
bal = mean(rec);

fprintf('A3 fine-tune complete: Acc=%.3f  BA=%.3f\n', acc, bal);
end
