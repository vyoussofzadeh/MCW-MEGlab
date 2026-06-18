%% ------------------------------------------------------------------------
% Step8_FineTune_Model_With_Feedback.m
%
% Fine-tune an existing trained model using reviewed feedback epochs.
% Keeps original validation/test set fixed for monitoring.
% -------------------------------------------------------------------------

clc; clear;

%% ======================== USER SETTINGS =================================

opts.modelFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Squiggles/Step5_lowmem_2DCNN_20260518_170934.mat';

opts.feedbackFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data/feedback_round01_epochs.mat';

opts.outModel = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Squiggles/Step8_model_feedback_round01_finetuned.mat';

opts.numEpochs = 5;
opts.miniBatchSize = 16;
opts.learnRate = 1e-5;      % small learning rate for fine-tuning

useGPU = canUseGPU;
fprintf('Use GPU: %d\n', useGPU);

rng(1);

%% ======================== LOAD MODEL ====================================

M = load(opts.modelFile);

if isfield(M, 'bestNet')
    net = M.bestNet;
elseif isfield(M, 'net')
    net = M.net;
else
    error('Model file does not contain bestNet or net.');
end

if isfield(M, 'datasetFile')
    datasetFile = M.datasetFile;
else
    error('Model file does not contain datasetFile.');
end

fprintf('Base model:\n%s\n', opts.modelFile);
fprintf('Base dataset:\n%s\n', datasetFile);

%% ======================== LOAD FEEDBACK =================================

F = load(opts.feedbackFile);

if isfield(F, 'feedback')
    Xfb = F.feedback.X;
    yfb = F.feedback.y(:);
else
    Xfb = F.X_feedback;
    yfb = F.y_feedback(:);
end

yfb = double(yfb(:));

fprintf('\nFeedback data:\n');
fprintf('Xfb: [%s]\n', num2str(size(Xfb)));
fprintf('True spikes:     %d\n', sum(yfb==1));
fprintf('False positives: %d\n', sum(yfb==0));

if numel(unique(yfb)) < 2
    warning('Feedback contains only one class. Fine-tuning may bias the model.');
end

C = size(Xfb,2);
T = size(Xfb,3);

%% ======================== LOAD ORIGINAL VAL/TEST =========================

[baseData, y_val, y_test] = open_base_dataset(datasetFile);

fprintf('\nOriginal validation/test:\n');
fprintf('Val labels:  spike=%d, nospike=%d\n', sum(y_val==1), sum(y_val==0));
fprintf('Test labels: spike=%d, nospike=%d\n', sum(y_test==1), sum(y_test==0));

%% ======================== BASELINE EVALUATION ============================

[valAcc0, valLoss0] = evaluateModelBase(net, baseData, 'X_val', y_val, opts.miniBatchSize, C, T, useGPU);
[testAcc0, testLoss0] = evaluateModelBase(net, baseData, 'X_test', y_test, opts.miniBatchSize, C, T, useGPU);

fprintf('\nBefore fine-tuning:\n');
fprintf('Val  acc %.2f%% | loss %.4f\n', 100*valAcc0, valLoss0);
fprintf('Test acc %.2f%% | loss %.4f\n', 100*testAcc0, testLoss0);

%% ======================== FINE-TUNING LOOP ===============================

bestNet = net;
bestValAcc = valAcc0;
bestValLoss = valLoss0;
bestEpoch = 0;

trailingAvg = [];
trailingAvgSq = [];
iteration = 0;

trainLossHistory = [];
valAccHistory = valAcc0;
valLossHistory = valLoss0;

Nfb = numel(yfb);

fprintf('\nStarting feedback fine-tuning...\n');

for epoch = 1:opts.numEpochs

    tic;

    order = randperm(Nfb);

    epochLoss = 0;
    nBatches = 0;

    for s = 1:opts.miniBatchSize:Nfb

        iteration = iteration + 1;

        e = min(s + opts.miniBatchSize - 1, Nfb);
        idx = order(s:e);

        [dlX, dlT] = make_dlX_dlT_for_model(net, Xfb(idx,:,:), yfb(idx), useGPU);

        [loss, gradients] = dlfeval(@modelGradients, net, dlX, dlT);

        [net, trailingAvg, trailingAvgSq] = adamupdate( ...
            net, gradients, trailingAvg, trailingAvgSq, iteration, opts.learnRate);

        epochLoss = epochLoss + double(gather(extractdata(loss)));
        nBatches = nBatches + 1;
    end

    meanTrainLoss = epochLoss / max(nBatches,1);

    [valAcc, valLoss] = evaluateModelBase(net, baseData, 'X_val', y_val, opts.miniBatchSize, C, T, useGPU);

    trainLossHistory(end+1,1) = meanTrainLoss;
    valAccHistory(end+1,1) = valAcc;
    valLossHistory(end+1,1) = valLoss;

    fprintf('Epoch %d/%d | feedback loss %.4f | val loss %.4f | val acc %.2f%% | %.1f sec\n', ...
        epoch, opts.numEpochs, meanTrainLoss, valLoss, 100*valAcc, toc);

    if valAcc > bestValAcc
        bestValAcc = valAcc;
        bestValLoss = valLoss;
        bestNet = net;
        bestEpoch = epoch;

        fprintf('*** New best fine-tuned model at epoch %d: val acc %.2f%% ***\n', ...
            epoch, 100*valAcc);
    end
end

%% ======================== FINAL EVALUATION ===============================

netFinal = net;
net = bestNet;

[valAcc, valLoss, yPredVal, pVal] = evaluateModelBase(net, baseData, 'X_val', y_val, opts.miniBatchSize, C, T, useGPU);
[testAcc, testLoss, yPredTest, pTest] = evaluateModelBase(net, baseData, 'X_test', y_test, opts.miniBatchSize, C, T, useGPU);

fprintf('\nBest fine-tuned epoch: %d\n', bestEpoch);
fprintf('Val  acc %.2f%% | loss %.4f\n', 100*valAcc, valLoss);
fprintf('Test acc %.2f%% | loss %.4f\n', 100*testAcc, testLoss);

%% ======================== SAVE UPDATED MODEL =============================

save(opts.outModel, ...
    'bestNet','netFinal','bestEpoch', ...
    'opts','datasetFile', ...
    'valAcc0','valLoss0','testAcc0','testLoss0', ...
    'valAcc','valLoss','testAcc','testLoss', ...
    'yPredVal','yPredTest','pVal','pTest', ...
    'trainLossHistory','valAccHistory','valLossHistory', ...
    '-v7.3');

fprintf('\nSaved fine-tuned model:\n%s\n', opts.outModel);

%% ========================================================================
% Helper functions
% ========================================================================

function [baseData, y_val, y_test] = open_base_dataset(datasetFile)

    info = whos('-file', datasetFile);
    names = {info.name};

    if all(ismember({'X_val','X_test','y_val','y_test'}, names))
        baseData.mode = 'matfile';
        baseData.m = matfile(datasetFile);

        S = load(datasetFile, 'y_val','y_test');
        y_val = double(S.y_val(:));
        y_test = double(S.y_test(:));

    elseif ismember('dataset', names)
        S = load(datasetFile, 'dataset');
        baseData.mode = 'struct';
        baseData.dataset = S.dataset;

        y_val = double(S.dataset.y_val(:));
        y_test = double(S.dataset.y_test(:));

    else
        error('Unknown dataset format: %s', datasetFile);
    end
end

function [dlX, dlT] = make_dlX_dlT_for_model(net, X, y, useGPU)

    % X: B x C x T

    inputLayer = net.Layers(1);
    inputClass = class(inputLayer);

    B = size(X,1);

    if contains(inputClass, 'ImageInputLayer')
        X = reshape(single(X), [B size(X,2) size(X,3) 1]);
        X = permute(X, [2 3 4 1]);      % C x T x 1 x B
        dlX = dlarray(X, 'SSCB');

    elseif contains(inputClass, 'SequenceInputLayer')
        X = permute(single(X), [2 3 1]); % C x T x B
        dlX = dlarray(X, 'CTB');

    else
        error('Unsupported input layer type: %s', inputClass);
    end

    Tmat = zeros(2, B, 'single');
    Tmat(1, y == 0) = 1;
    Tmat(2, y == 1) = 1;

    dlT = dlarray(Tmat, 'CB');

    if useGPU
        dlX = gpuArray(dlX);
        dlT = gpuArray(dlT);
    end
end

function [loss, gradients] = modelGradients(net, dlX, dlT)

    dlY = forward(net, dlX);

    epsVal = 1e-7;
    loss = -mean(sum(dlT .* log(dlY + epsVal), 1));

    gradients = dlgradient(loss, net.Learnables);
end

function [acc, meanLoss, yPredAll, pSpikeAll] = evaluateModelBase(net, baseData, xVar, y, miniBatchSize, C, T, useGPU)

    N = numel(y);

    yPredAll = zeros(N,1);
    pSpikeAll = zeros(N,1);

    lossSum = 0;
    nBatches = 0;

    for s = 1:miniBatchSize:N

        e = min(s + miniBatchSize - 1, N);
        idx = s:e;

        X = read_base_X(baseData, xVar, idx, C, T);

        [dlX, dlT] = make_dlX_dlT_for_model(net, X, y(idx), useGPU);

        dlY = predict(net, dlX);

        epsVal = 1e-7;
        loss = -mean(sum(dlT .* log(dlY + epsVal), 1));

        P = gather(extractdata(dlY));

        [~, pc] = max(P, [], 1);
        yPred = pc(:) - 1;

        yPredAll(idx) = yPred;
        pSpikeAll(idx) = P(2,:)';

        lossSum = lossSum + double(gather(extractdata(loss)));
        nBatches = nBatches + 1;
    end

    acc = mean(yPredAll == y);
    meanLoss = lossSum / max(nBatches,1);
end

function X = read_base_X(baseData, xVar, idx, C, T)

    switch baseData.mode

        case 'matfile'
            switch xVar
                case 'X_val'
                    X = baseData.m.X_val(idx,:,:);
                case 'X_test'
                    X = baseData.m.X_test(idx,:,:);
                otherwise
                    error('Unsupported xVar: %s', xVar);
            end

        case 'struct'
            switch xVar
                case 'X_val'
                    X = baseData.dataset.X_val(idx,:,:);
                case 'X_test'
                    X = baseData.dataset.X_test(idx,:,:);
                otherwise
                    error('Unsupported xVar: %s', xVar);
            end
    end

    if size(X,1) ~= numel(idx)
        X = reshape(X, [numel(idx) C T]);
    end

    X = single(X);
end