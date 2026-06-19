%% ------------------------------------------------------------------------
% Step8_FineTune_Model_With_Feedback_FAST.m
%
% CPU-friendly fine-tuning using Step7b feedback epochs.
%
% Faster than the full Step8 because it:
%   - skips baseline test by default
%   - uses larger batches for evaluation
%   - tests only once at the end
%   - optionally evaluates validation only every N epochs
%
% Input:
%   Step7b feedback_round01_epochs.mat
%
% Output:
%   fine-tuned model MAT with bestNet and netFinal
% -------------------------------------------------------------------------

clc; clear;

%% ======================== USER SETTINGS =================================

opts.feedbackFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data/feedback_round01_epochs.mat';

% Leave empty to use feedback.modelFile from Step7b.
opts.modelFile = '';

% Leave empty to use datasetFile saved in model/feedback.
opts.datasetFileOverride = '';

opts.outModel = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Squiggles/Step8_model_feedback_round01_finetuned_FAST.mat';

% Fine-tuning settings. Keep this short for 74 reviewed epochs.
opts.numEpochs = 3;
opts.trainMiniBatchSize = 16;
opts.evalMiniBatchSize  = 256;
opts.learnRate = 1e-5;

% Evaluation settings.
opts.doBaselineFeedback = true;
opts.doBaselineVal      = true;
opts.doBaselineTest     = false;
opts.doFinalFeedback    = true;
opts.doFinalVal         = true;
opts.doFinalTest        = true;

% If CPU is slow, set to 2 or 3. Validation is still run on final best model.
opts.validateEveryEpoch = 1;

% Best model selection.
opts.bestMetric = 'valloss';       % 'valloss' or 'valacc'

% Your older Step8 assumed the network already returns softmax probabilities.
% Leave false unless your dlnetwork returns logits.
opts.applySoftmaxToOutput = false;

useGPU = canUseGPU;
fprintf('Use GPU: %d\n', useGPU);

rng(1);

%% ======================== LOAD FEEDBACK =================================

fprintf('\nLoading feedback:\n%s\n', opts.feedbackFile);

F = load(opts.feedbackFile);

if isfield(F, 'feedback')
    feedback = F.feedback;
    Xfb = feedback.X;
    yfb = feedback.y(:);
else
    feedback = struct();
    Xfb = F.X_feedback;
    yfb = F.y_feedback(:);
end

Xfb = single(Xfb);
yfb = double(yfb(:));

keep = yfb == 0 | yfb == 1;
Xfb = Xfb(keep,:,:);
yfb = yfb(keep);

if isempty(yfb)
    error('Feedback file has no usable labels 0/1.');
end

fprintf('\nFeedback data:\n');
fprintf('Xfb: [%s]\n', num2str(size(Xfb)));
fprintf('True spikes:     %d\n', sum(yfb==1));
fprintf('False positives: %d\n', sum(yfb==0));

if numel(unique(yfb)) < 2
    warning('Feedback contains only one class. Fine-tuning may bias the model.');
end

C = size(Xfb,2);
T = size(Xfb,3);

%% ======================== LOAD MODEL ====================================

if isempty(opts.modelFile)
    if isfield(feedback, 'modelFile') && ~isempty(feedback.modelFile)
        opts.modelFile = feedback.modelFile;
    else
        error('Set opts.modelFile, or use a Step7b file with feedback.modelFile.');
    end
end

fprintf('\nLoading base model:\n%s\n', opts.modelFile);

M = load(opts.modelFile);

if isfield(M, 'bestNet')
    net = M.bestNet;
elseif isfield(M, 'net')
    net = M.net;
else
    error('Model file does not contain bestNet or net.');
end

baseNet = net;

datasetFile = '';
if ~isempty(opts.datasetFileOverride)
    datasetFile = opts.datasetFileOverride;
elseif isfield(M, 'datasetFile') && ~isempty(M.datasetFile)
    datasetFile = M.datasetFile;
elseif isfield(feedback, 'datasetFile') && ~isempty(feedback.datasetFile)
    datasetFile = feedback.datasetFile;
end

datasetFile = char(datasetFile);

if isempty(datasetFile) || ~isfile(datasetFile)
    error('Could not find original datasetFile. Set opts.datasetFileOverride.');
end

fprintf('Base dataset:\n%s\n', datasetFile);

%% ======================== LOAD ORIGINAL VAL/TEST =========================

[baseData, y_val, y_test, szVal, szTest] = open_base_dataset(datasetFile);

fprintf('\nOriginal validation/test:\n');
fprintf('X_val:  [%s]\n', num2str(szVal));
fprintf('X_test: [%s]\n', num2str(szTest));
fprintf('Val labels:  spike=%d, nospike=%d\n', sum(y_val==1), sum(y_val==0));
fprintf('Test labels: spike=%d, nospike=%d\n', sum(y_test==1), sum(y_test==0));

if szVal(2) ~= C || szVal(3) ~= T
    error('Feedback size [%d x %d] does not match validation size [%d x %d].', ...
        C, T, szVal(2), szVal(3));
end

%% ======================== BASELINE EVALUATION ============================

fbAcc0 = NaN; fbLoss0 = NaN; yPredFb0 = []; pFb0 = [];
valAcc0 = NaN; valLoss0 = NaN;
testAcc0 = NaN; testLoss0 = NaN;

if opts.doBaselineFeedback
    fprintf('\nEvaluating base model on feedback set...\n');
    [fbAcc0, fbLoss0, yPredFb0, pFb0] = evaluateArray( ...
        baseNet, Xfb, yfb, opts.evalMiniBatchSize, useGPU, opts.applySoftmaxToOutput);
end

if opts.doBaselineVal
    fprintf('\nEvaluating base model on validation set...\n');
    [valAcc0, valLoss0] = evaluateBase( ...
        baseNet, baseData, 'X_val', y_val, opts.evalMiniBatchSize, C, T, useGPU, opts.applySoftmaxToOutput);
end

if opts.doBaselineTest
    fprintf('\nEvaluating base model on test set...\n');
    [testAcc0, testLoss0] = evaluateBase( ...
        baseNet, baseData, 'X_test', y_test, opts.evalMiniBatchSize, C, T, useGPU, opts.applySoftmaxToOutput);
end

fprintf('\nBefore fine-tuning:\n');
fprintf('Feedback acc %.2f%% | loss %.4f\n', 100*fbAcc0, fbLoss0);
fprintf('Val      acc %.2f%% | loss %.4f\n', 100*valAcc0, valLoss0);
if opts.doBaselineTest
    fprintf('Test     acc %.2f%% | loss %.4f\n', 100*testAcc0, testLoss0);
else
    fprintf('Test baseline skipped.\n');
end

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
fprintf('numEpochs=%d | trainBatch=%d | evalBatch=%d | learnRate=%g\n', ...
    opts.numEpochs, opts.trainMiniBatchSize, opts.evalMiniBatchSize, opts.learnRate);

for epoch = 1:opts.numEpochs

    tic;

    order = randperm(Nfb);

    epochLoss = 0;
    nBatches = 0;

    for s = 1:opts.trainMiniBatchSize:Nfb

        iteration = iteration + 1;

        e = min(s + opts.trainMiniBatchSize - 1, Nfb);
        idx = order(s:e);

        [dlX, dlT] = make_dlX_dlT_for_model(net, Xfb(idx,:,:), yfb(idx), useGPU);

        [loss, gradients] = dlfeval( ...
            @modelGradients, net, dlX, dlT, opts.applySoftmaxToOutput);

        [net, trailingAvg, trailingAvgSq] = adamupdate( ...
            net, gradients, trailingAvg, trailingAvgSq, iteration, opts.learnRate);

        epochLoss = epochLoss + double(gather(extractdata(loss)));
        nBatches = nBatches + 1;
    end

    meanTrainLoss = epochLoss / max(nBatches,1);

    doValidateNow = mod(epoch, opts.validateEveryEpoch) == 0 || epoch == opts.numEpochs;

    if doValidateNow
        fprintf('Evaluating validation after epoch %d...\n', epoch);
        [valAcc, valLoss] = evaluateBase( ...
            net, baseData, 'X_val', y_val, opts.evalMiniBatchSize, C, T, useGPU, opts.applySoftmaxToOutput);
    else
        valAcc = NaN;
        valLoss = NaN;
    end

    trainLossHistory(end+1,1) = meanTrainLoss;
    valAccHistory(end+1,1) = valAcc;
    valLossHistory(end+1,1) = valLoss;

    fprintf('Epoch %d/%d | feedback loss %.4f | val loss %.4f | val acc %.2f%% | %.1f sec\n', ...
        epoch, opts.numEpochs, meanTrainLoss, valLoss, 100*valAcc, toc);

    if doValidateNow && is_better_model(valAcc, valLoss, bestValAcc, bestValLoss, opts.bestMetric)
        bestValAcc = valAcc;
        bestValLoss = valLoss;
        bestNet = net;
        bestEpoch = epoch;

        fprintf('*** New best model at epoch %d: val acc %.2f%% | val loss %.4f ***\n', ...
            epoch, 100*valAcc, valLoss);
    end
end

%% ======================== FINAL EVALUATION ===============================

netFinal = net;
net = bestNet;

fbAcc = NaN; fbLoss = NaN; yPredFb = []; pFb = [];
valAcc = NaN; valLoss = NaN; yPredVal = []; pVal = [];
testAcc = NaN; testLoss = NaN; yPredTest = []; pTest = [];

if opts.doFinalFeedback
    fprintf('\nFinal evaluation on feedback set...\n');
    [fbAcc, fbLoss, yPredFb, pFb] = evaluateArray( ...
        net, Xfb, yfb, opts.evalMiniBatchSize, useGPU, opts.applySoftmaxToOutput);
end

if opts.doFinalVal
    fprintf('\nFinal evaluation on validation set...\n');
    [valAcc, valLoss, yPredVal, pVal] = evaluateBase( ...
        net, baseData, 'X_val', y_val, opts.evalMiniBatchSize, C, T, useGPU, opts.applySoftmaxToOutput);
end

if opts.doFinalTest
    fprintf('\nFinal evaluation on test set...\n');
    [testAcc, testLoss, yPredTest, pTest] = evaluateBase( ...
        net, baseData, 'X_test', y_test, opts.evalMiniBatchSize, C, T, useGPU, opts.applySoftmaxToOutput);
end

fprintf('\nBest fine-tuned epoch: %d\n', bestEpoch);
fprintf('Feedback acc %.2f%% -> %.2f%% | loss %.4f -> %.4f\n', ...
    100*fbAcc0, 100*fbAcc, fbLoss0, fbLoss);
fprintf('Val      acc %.2f%% -> %.2f%% | loss %.4f -> %.4f\n', ...
    100*valAcc0, 100*valAcc, valLoss0, valLoss);
fprintf('Test     acc %.2f%% | loss %.4f\n', 100*testAcc, testLoss);

%% ======================== SAVE UPDATED MODEL =============================

save(opts.outModel, ...
    'bestNet','netFinal','baseNet','bestEpoch', ...
    'opts','feedback','datasetFile', ...
    'fbAcc0','fbLoss0','valAcc0','valLoss0','testAcc0','testLoss0', ...
    'fbAcc','fbLoss','valAcc','valLoss','testAcc','testLoss', ...
    'yPredFb0','pFb0','yPredFb','pFb', ...
    'yPredVal','yPredTest','pVal','pTest', ...
    'trainLossHistory','valAccHistory','valLossHistory', ...
    '-v7.3');

fprintf('\nSaved fine-tuned model:\n%s\n', opts.outModel);

%% ========================================================================
% Helper functions
% ========================================================================

function [baseData, y_val, y_test, szVal, szTest] = open_base_dataset(datasetFile)

    info = whos('-file', datasetFile);
    names = {info.name};

    if all(ismember({'X_val','X_test','y_val','y_test'}, names))

        baseData.mode = 'matfile';
        baseData.m = matfile(datasetFile);

        szVal = size(baseData.m, 'X_val');
        szTest = size(baseData.m, 'X_test');

        S = load(datasetFile, 'y_val','y_test');
        y_val = double(S.y_val(:));
        y_test = double(S.y_test(:));

    elseif ismember('dataset', names)

        S = load(datasetFile, 'dataset');

        baseData.mode = 'struct';
        baseData.dataset = S.dataset;

        szVal = size(S.dataset.X_val);
        szTest = size(S.dataset.X_test);

        y_val = double(S.dataset.y_val(:));
        y_test = double(S.dataset.y_test(:));

    else
        error('Unknown dataset format: %s', datasetFile);
    end
end

function [dlX, dlT] = make_dlX_dlT_for_model(net, X, y, useGPU)

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

function [loss, gradients] = modelGradients(net, dlX, dlT, applySoftmaxToOutput)

    dlY = forward(net, dlX);

    if applySoftmaxToOutput
        dlY = softmax(dlY);
    end

    epsVal = 1e-7;
    loss = -mean(sum(dlT .* log(dlY + epsVal), 1));

    gradients = dlgradient(loss, net.Learnables);
end

function [acc, meanLoss, yPredAll, pSpikeAll] = evaluateBase( ...
    net, baseData, xVar, y, miniBatchSize, C, T, useGPU, applySoftmaxToOutput)

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

        if applySoftmaxToOutput
            dlY = softmax(dlY);
        end

        epsVal = 1e-7;
        loss = -mean(sum(dlT .* log(dlY + epsVal), 1));

        P = gather(extractdata(dlY));

        [~, pc] = max(P, [], 1);
        yPred = pc(:) - 1;

        yPredAll(idx) = yPred;
        pSpikeAll(idx) = P(2,:)';

        lossSum = lossSum + double(gather(extractdata(loss)));
        nBatches = nBatches + 1;

        if mod(nBatches, 10) == 0
            fprintf('  eval %s: %d / %d\n', xVar, e, N);
        end
    end

    acc = mean(yPredAll == y);
    meanLoss = lossSum / max(nBatches,1);
end

function [acc, meanLoss, yPredAll, pSpikeAll] = evaluateArray( ...
    net, Xall, y, miniBatchSize, useGPU, applySoftmaxToOutput)

    N = numel(y);

    yPredAll = zeros(N,1);
    pSpikeAll = zeros(N,1);

    lossSum = 0;
    nBatches = 0;

    for s = 1:miniBatchSize:N

        e = min(s + miniBatchSize - 1, N);
        idx = s:e;

        [dlX, dlT] = make_dlX_dlT_for_model(net, Xall(idx,:,:), y(idx), useGPU);

        dlY = predict(net, dlX);

        if applySoftmaxToOutput
            dlY = softmax(dlY);
        end

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

        otherwise
            error('Unknown baseData mode: %s', baseData.mode);
    end

    if size(X,1) ~= numel(idx)
        X = reshape(X, [numel(idx) C T]);
    end

    X = single(X);
end

function tf = is_better_model(valAcc, valLoss, bestValAcc, bestValLoss, metric)

    if isnan(bestValAcc) || isnan(bestValLoss)
        tf = true;
        return;
    end

    switch lower(metric)
        case 'valacc'
            tf = valAcc > bestValAcc || ...
                (abs(valAcc - bestValAcc) < 1e-12 && valLoss < bestValLoss);
        case 'valloss'
            tf = valLoss < bestValLoss || ...
                (abs(valLoss - bestValLoss) < 1e-12 && valAcc > bestValAcc);
        otherwise
            error('Unknown opts.bestMetric: %s', metric);
    end
end