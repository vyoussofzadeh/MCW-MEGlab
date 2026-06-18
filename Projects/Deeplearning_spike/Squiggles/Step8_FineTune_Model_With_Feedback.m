%% ------------------------------------------------------------------------
% Step8_FineTune_Model_With_Feedback.m
%
% Fine-tune an existing trained model using Step7b reviewed feedback epochs.
% Keeps the original validation/test sets fixed for monitoring.
%
% Recommended flow:
%   Step7a -> *_interactive_reviewed.mat
%   Step7b -> feedback_round01_epochs.mat
%   Step8  -> fine-tuned model MAT
%
% This version can mix reviewed feedback with a small random replay sample
% from the original training set to reduce catastrophic forgetting.
% -------------------------------------------------------------------------

clc; clear;

%% ======================== USER SETTINGS =================================

% Step7b output.
opts.feedbackFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data/feedback_round01_epochs.mat';

% Leave empty to use feedback.modelFile saved by Step7b.
opts.modelFile = '';

% Leave empty to use datasetFile saved in the model / feedback file.
opts.datasetFileOverride = '';

opts.outModel = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Squiggles/Step8_model_feedback_round01_finetuned.mat';

opts.numEpochs = 8;
opts.miniBatchSize = 16;
opts.learnRate = 1e-5;

% Replay original training examples during fine-tuning.
% This is usually safer than training only on reviewed candidates.
opts.useReplay = true;
opts.replayPerFeedback = 1.0;       % 1.0 means one original train epoch per feedback epoch size
opts.maxReplayPerEpoch = 512;       % cap for speed/memory
opts.balanceReplayClasses = true;

% If true, duplicate the minority feedback class so each feedback epoch is
% approximately balanced before adding replay examples.
opts.balanceFeedbackClasses = true;

% Best model selection. 'valloss' is usually steadier than 'valacc'.
opts.bestMetric = 'valloss';        % 'valloss' or 'valacc'

% Your older Step8 assumed the network already ends with softmax and returns
% probabilities. Leave this false for that case. Set true only if your
% dlnetwork returns raw logits.
opts.applySoftmaxToOutput = false;

useGPU = canUseGPU;
fprintf('Use GPU: %d\n', useGPU);

rng(1);

%% ======================== LOAD FEEDBACK =================================

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

yfb = double(yfb(:));

keep = yfb == 0 | yfb == 1;
Xfb = Xfb(keep,:,:);
yfb = yfb(keep);

if isempty(yfb)
    error('Feedback file has no usable labels 0/1: %s', opts.feedbackFile);
end

fprintf('\nFeedback data:\n');
fprintf('File: %s\n', opts.feedbackFile);
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
        error('Set opts.modelFile, or provide a Step7b feedback file with feedback.modelFile.');
    end
end

M = load(opts.modelFile);

if isfield(M, 'bestNet')
    net = M.bestNet;
elseif isfield(M, 'net')
    net = M.net;
else
    error('Model file does not contain bestNet or net.');
end

if ~isa(net, 'dlnetwork')
    error(['This Step8 script expects a dlnetwork saved as bestNet/net.\n' ...
           'Loaded model class is: %s'], class(net));
end

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

fprintf('\nBase model:\n%s\n', opts.modelFile);
fprintf('\nBase dataset:\n%s\n', datasetFile);

%% ======================== LOAD ORIGINAL DATASET ==========================

[baseData, y_train, y_val, y_test, Cbase, Tbase] = open_base_dataset(datasetFile);

if Cbase ~= C || Tbase ~= T
    error('Feedback size is C=%d,T=%d but base dataset is C=%d,T=%d.', C, T, Cbase, Tbase);
end

fprintf('\nOriginal dataset labels:\n');
fprintf('Train labels: spike=%d, nospike=%d\n', sum(y_train==1), sum(y_train==0));
fprintf('Val labels:   spike=%d, nospike=%d\n', sum(y_val==1), sum(y_val==0));
fprintf('Test labels:  spike=%d, nospike=%d\n', sum(y_test==1), sum(y_test==0));

%% ======================== BASELINE EVALUATION ============================

[valAcc0, valLoss0, yPredVal0, pVal0] = evaluateModelBase( ...
    net, baseData, 'X_val', y_val, opts.miniBatchSize, C, T, useGPU);

[testAcc0, testLoss0, yPredTest0, pTest0] = evaluateModelBase( ...
    net, baseData, 'X_test', y_test, opts.miniBatchSize, C, T, useGPU);

[fbAcc0, fbLoss0, yPredFb0, pFb0] = evaluateModelArray( ...
    net, Xfb, yfb, opts.miniBatchSize, useGPU);

fprintf('\nBefore fine-tuning:\n');
fprintf('Feedback acc %.2f%% | loss %.4f\n', 100*fbAcc0, fbLoss0);
fprintf('Val      acc %.2f%% | loss %.4f\n', 100*valAcc0, valLoss0);
fprintf('Test     acc %.2f%% | loss %.4f\n', 100*testAcc0, testLoss0);

%% ======================== FINE-TUNING LOOP ===============================

bestNet = net;
bestValAcc = valAcc0;
bestValLoss = valLoss0;
bestEpoch = 0;

trailingAvg = [];
trailingAvgSq = [];
iteration = 0;

trainLossHistory = [];
feedbackAccHistory = fbAcc0;
feedbackLossHistory = fbLoss0;
valAccHistory = valAcc0;
valLossHistory = valLoss0;
testAccHistory = testAcc0;
testLossHistory = testLoss0;

fprintf('\nStarting feedback fine-tuning...\n');

for epoch = 1:opts.numEpochs

    tic;

    idxFb = make_feedback_epoch_indices(yfb, opts.balanceFeedbackClasses);

    Xepoch = Xfb(idxFb,:,:);
    yepoch = yfb(idxFb);

    if opts.useReplay
        nReplay = min(round(numel(idxFb) * opts.replayPerFeedback), opts.maxReplayPerEpoch);
        idxReplay = make_replay_indices(y_train, nReplay, opts.balanceReplayClasses);

        Xreplay = read_base_X(baseData, 'X_train', idxReplay, C, T);
        yreplay = y_train(idxReplay);

        Xepoch = cat(1, Xepoch, Xreplay);
        yepoch = [yepoch; yreplay(:)];
    else
        idxReplay = [];
    end

    Nepoch = numel(yepoch);
    order = randperm(Nepoch);

    epochLoss = 0;
    nBatches = 0;

    for s = 1:opts.miniBatchSize:Nepoch

        iteration = iteration + 1;

        e = min(s + opts.miniBatchSize - 1, Nepoch);
        idx = order(s:e);

        [dlX, dlT] = make_dlX_dlT_for_model(net, Xepoch(idx,:,:), yepoch(idx), useGPU);

        [loss, gradients] = dlfeval(@modelGradients, net, dlX, dlT, opts.applySoftmaxToOutput);

        [net, trailingAvg, trailingAvgSq] = adamupdate( ...
            net, gradients, trailingAvg, trailingAvgSq, iteration, opts.learnRate);

        epochLoss = epochLoss + double(gather(extractdata(loss)));
        nBatches = nBatches + 1;
    end

    meanTrainLoss = epochLoss / max(nBatches,1);

    [fbAcc, fbLoss] = evaluateModelArray(net, Xfb, yfb, opts.miniBatchSize, useGPU);
    [valAcc, valLoss] = evaluateModelBase(net, baseData, 'X_val', y_val, opts.miniBatchSize, C, T, useGPU);
    [testAcc, testLoss] = evaluateModelBase(net, baseData, 'X_test', y_test, opts.miniBatchSize, C, T, useGPU);

    trainLossHistory(end+1,1) = meanTrainLoss;
    feedbackAccHistory(end+1,1) = fbAcc;
    feedbackLossHistory(end+1,1) = fbLoss;
    valAccHistory(end+1,1) = valAcc;
    valLossHistory(end+1,1) = valLoss;
    testAccHistory(end+1,1) = testAcc;
    testLossHistory(end+1,1) = testLoss;

    fprintf(['Epoch %d/%d | train %.4f | feedback %.2f%%/%.4f | ' ...
             'val %.2f%%/%.4f | test %.2f%%/%.4f | replay=%d | %.1f sec\n'], ...
        epoch, opts.numEpochs, meanTrainLoss, ...
        100*fbAcc, fbLoss, 100*valAcc, valLoss, 100*testAcc, testLoss, ...
        numel(idxReplay), toc);

    if is_better_model(valAcc, valLoss, bestValAcc, bestValLoss, opts.bestMetric)
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

[fbAcc, fbLoss, yPredFb, pFb] = evaluateModelArray( ...
    net, Xfb, yfb, opts.miniBatchSize, useGPU);

[valAcc, valLoss, yPredVal, pVal] = evaluateModelBase( ...
    net, baseData, 'X_val', y_val, opts.miniBatchSize, C, T, useGPU);

[testAcc, testLoss, yPredTest, pTest] = evaluateModelBase( ...
    net, baseData, 'X_test', y_test, opts.miniBatchSize, C, T, useGPU);

fprintf('\nBest fine-tuned epoch: %d\n', bestEpoch);
fprintf('Feedback acc %.2f%% | loss %.4f\n', 100*fbAcc, fbLoss);
fprintf('Val      acc %.2f%% | loss %.4f\n', 100*valAcc, valLoss);
fprintf('Test     acc %.2f%% | loss %.4f\n', 100*testAcc, testLoss);

%% ======================== SAVE UPDATED MODEL =============================

save(opts.outModel, ...
    'bestNet','netFinal','bestEpoch', ...
    'opts','feedback','datasetFile', ...
    'valAcc0','valLoss0','testAcc0','testLoss0','fbAcc0','fbLoss0', ...
    'valAcc','valLoss','testAcc','testLoss','fbAcc','fbLoss', ...
    'yPredVal0','yPredTest0','yPredFb0','pVal0','pTest0','pFb0', ...
    'yPredVal','yPredTest','yPredFb','pVal','pTest','pFb', ...
    'trainLossHistory','feedbackAccHistory','feedbackLossHistory', ...
    'valAccHistory','valLossHistory','testAccHistory','testLossHistory', ...
    '-v7.3');

fprintf('\nSaved fine-tuned model:\n%s\n', opts.outModel);

%% ========================================================================
% Helper functions
% ========================================================================

function [baseData, y_train, y_val, y_test, C, T] = open_base_dataset(datasetFile)

    info = whos('-file', datasetFile);
    names = {info.name};

    if all(ismember({'X_train','X_val','X_test','y_train','y_val','y_test'}, names))
        baseData.mode = 'matfile';
        baseData.m = matfile(datasetFile);

        sz = size(baseData.m, 'X_train');
        C = sz(2);
        T = sz(3);

        S = load(datasetFile, 'y_train','y_val','y_test');
        y_train = double(S.y_train(:));
        y_val = double(S.y_val(:));
        y_test = double(S.y_test(:));

    elseif ismember('dataset', names)
        S = load(datasetFile, 'dataset');
        baseData.mode = 'struct';
        baseData.dataset = S.dataset;

        sz = size(S.dataset.X_train);
        C = sz(2);
        T = sz(3);

        y_train = double(S.dataset.y_train(:));
        y_val = double(S.dataset.y_val(:));
        y_test = double(S.dataset.y_test(:));

    else
        error('Unknown dataset format: %s', datasetFile);
    end
end

function idx = make_feedback_epoch_indices(y, doBalance)

    y = double(y(:));

    if ~doBalance || numel(unique(y)) < 2
        idx = (1:numel(y))';
        return;
    end

    idx0 = find(y == 0);
    idx1 = find(y == 1);
    n = max(numel(idx0), numel(idx1));

    idx0 = idx0(randi(numel(idx0), n, 1));
    idx1 = idx1(randi(numel(idx1), n, 1));

    idx = [idx0; idx1];
end

function idx = make_replay_indices(y_train, nReplay, doBalance)

    y_train = double(y_train(:));

    if nReplay <= 0
        idx = [];
        return;
    end

    if ~doBalance || numel(unique(y_train)) < 2
        idx = randi(numel(y_train), nReplay, 1);
        return;
    end

    n0 = floor(nReplay/2);
    n1 = nReplay - n0;

    idx0all = find(y_train == 0);
    idx1all = find(y_train == 1);

    idx0 = idx0all(randi(numel(idx0all), n0, 1));
    idx1 = idx1all(randi(numel(idx1all), n1, 1));

    idx = [idx0; idx1];
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

function [acc, meanLoss, yPredAll, pSpikeAll] = evaluateModelArray(net, Xall, y, miniBatchSize, useGPU)

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

        if opts_apply_softmax_from_base_workspace()
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

        if opts_apply_softmax_from_base_workspace()
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

function tf = opts_apply_softmax_from_base_workspace()

    % Nested-free helper functions cannot see opts directly. This keeps the
    % evaluation helpers in sync with the script setting without changing all
    % historical call signatures.
    try
        optsLocal = evalin('caller', 'opts');
        tf = isfield(optsLocal, 'applySoftmaxToOutput') && optsLocal.applySoftmaxToOutput;
    catch
        tf = false;
    end
end

function tf = is_better_model(valAcc, valLoss, bestValAcc, bestValLoss, metric)

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

function X = read_base_X(baseData, xVar, idx, C, T)

    switch baseData.mode

        case 'matfile'
            switch xVar
                case 'X_train'
                    X = baseData.m.X_train(idx,:,:);
                case 'X_val'
                    X = baseData.m.X_val(idx,:,:);
                case 'X_test'
                    X = baseData.m.X_test(idx,:,:);
                otherwise
                    error('Unsupported xVar: %s', xVar);
            end

        case 'struct'
            switch xVar
                case 'X_train'
                    X = baseData.dataset.X_train(idx,:,:);
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