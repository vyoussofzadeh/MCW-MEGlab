%% ------------------------------------------------------------------------
% Step5_QC_and_train_baseline_model_lowmem.m
%
% Low-memory training from dataset_DL_ready.mat.
% Does NOT load X_train/X_val/X_test fully into RAM.
% Reads mini-batches from matfile.
% -------------------------------------------------------------------------

clc; clear;

% datasetFile = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files/dataset_DL_ready.mat';
datasetFile = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files/dataset_DL_ready_win_m250_p500.mat';

%% ---- Inspect saved variables -------------------------------------------

whos('-file', datasetFile)

m = matfile(datasetFile);

szTrain = size(m, 'X_train');
szVal   = size(m, 'X_val');
szTest  = size(m, 'X_test');

fprintf('\nArray sizes:\n');
fprintf('X_train: [%s]\n', num2str(szTrain));
fprintf('X_val:   [%s]\n', num2str(szVal));
fprintf('X_test:  [%s]\n', num2str(szTest));

S = load(datasetFile, ...
    'dataset', ...
    'y_train','y_val','y_test', ...
    'subject_id_train','subject_id_val','subject_id_test', ...
    'base_train','base_val','base_test');

dataset = S.dataset;

y_train = S.y_train(:);
y_val   = S.y_val(:);
y_test  = S.y_test(:);

fprintf('\nClass counts:\n');
fprintf('Train: spike=%d, nospike=%d\n', sum(y_train==1), sum(y_train==0));
fprintf('Val:   spike=%d, nospike=%d\n', sum(y_val==1),   sum(y_val==0));
fprintf('Test:  spike=%d, nospike=%d\n', sum(y_test==1),  sum(y_test==0));

fprintf('\nSubject leakage check:\n');

trainSubj = unique(S.subject_id_train);
valSubj   = unique(S.subject_id_val);
testSubj  = unique(S.subject_id_test);

fprintf('Train-Val overlap:  %d\n', numel(intersect(trainSubj, valSubj)));
fprintf('Train-Test overlap: %d\n', numel(intersect(trainSubj, testSubj)));
fprintf('Val-Test overlap:   %d\n', numel(intersect(valSubj, testSubj)));

%% ---- Settings -----------------------------------------------------------

Ntrain = szTrain(1);
C      = szTrain(2);
T      = szTrain(3);

numClasses = 2;

miniBatchSize = 32;     % start safe; later try 64
numEpochs     = 15;
learnRate     = 1e-3;

useGPU = canUseGPU;
fprintf('\nUse GPU: %d\n', useGPU);

rng(1);

%% ---- Define a simple low-memory 2D CNN baseline -------------------------
% Input is treated as a 306 x 306 image:
%   rows = MEG channels
%   columns = time samples
%
% This is a baseline. Later we can switch to CNN1D/LSTM/datastore.

layers = [
    imageInputLayer([C T 1], 'Normalization','none', 'Name','input')

    convolution2dLayer([1 7], 16, 'Padding','same', 'Name','conv_time1')
    batchNormalizationLayer('Name','bn1')
    reluLayer('Name','relu1')
    maxPooling2dLayer([1 2], 'Stride',[1 2], 'Name','pool_time1')

    convolution2dLayer([3 5], 32, 'Padding','same', 'Name','conv_space_time2')
    batchNormalizationLayer('Name','bn2')
    reluLayer('Name','relu2')
    maxPooling2dLayer([2 2], 'Stride',[2 2], 'Name','pool2')

    convolution2dLayer([3 3], 64, 'Padding','same', 'Name','conv3')
    batchNormalizationLayer('Name','bn3')
    reluLayer('Name','relu3')

    globalAveragePooling2dLayer('Name','gap')

    fullyConnectedLayer(numClasses, 'Name','fc')
    softmaxLayer('Name','softmax')
];

net = dlnetwork(layers);

%% ---- Training loop ------------------------------------------------------

trailingAvg = [];
trailingAvgSq = [];
iteration = 0;

trainLossHistory = [];
valAccHistory = [];

fprintf('\nStarting low-memory training...\n');

for epoch = 1:numEpochs

    tic;

    order = randperm(Ntrain);

    epochLoss = 0;
    nBatches = 0;

    for bStart = 1:miniBatchSize:Ntrain

        iteration = iteration + 1;

        bEnd = min(bStart + miniBatchSize - 1, Ntrain);
        batchIdx = order(bStart:bEnd);

        [dlX, dlT] = readBatchFromMatfile(m, 'X_train', y_train, batchIdx, C, T, useGPU);

        [loss, gradients, state] = dlfeval(@modelGradients, net, dlX, dlT);
        net.State = state;

        [net, trailingAvg, trailingAvgSq] = adamupdate( ...
            net, gradients, trailingAvg, trailingAvgSq, iteration, learnRate);

        epochLoss = epochLoss + double(gather(extractdata(loss)));
        nBatches = nBatches + 1;

        if mod(nBatches, 25) == 0
            fprintf('Epoch %d/%d | batch %d | loss %.4f\n', ...
                epoch, numEpochs, nBatches, double(gather(extractdata(loss))));
        end
    end

    meanTrainLoss = epochLoss / max(nBatches,1);

    [valAcc, valLoss] = evaluateModelLowMem(net, m, 'X_val', y_val, miniBatchSize, C, T, useGPU);

    trainLossHistory(end+1,1) = meanTrainLoss; %#ok<SAGROW>
    valAccHistory(end+1,1) = valAcc; %#ok<SAGROW>

    fprintf('\nEpoch %d/%d complete | train loss %.4f | val loss %.4f | val acc %.2f%% | %.1f sec\n\n', ...
        epoch, numEpochs, meanTrainLoss, valLoss, 100*valAcc, toc);
end

%% ---- Final evaluation ---------------------------------------------------

[valAcc, valLoss, yPredVal] = evaluateModelLowMem(net, m, 'X_val', y_val, miniBatchSize, C, T, useGPU);
[testAcc, testLoss, yPredTest] = evaluateModelLowMem(net, m, 'X_test', y_test, miniBatchSize, C, T, useGPU);

fprintf('\nFinal validation accuracy: %.2f%% | loss %.4f\n', 100*valAcc, valLoss);
fprintf('Final test accuracy:       %.2f%% | loss %.4f\n', 100*testAcc, testLoss);

figure;
confusionchart(categorical(y_test, [0 1], {'NoSpike','Spike'}), ...
               categorical(yPredTest, [0 1], {'NoSpike','Spike'}));
title('Test confusion matrix');

figure;
plot(trainLossHistory, '-o');
xlabel('Epoch');
ylabel('Train loss');
title('Training loss');
grid on;

figure;
plot(100*valAccHistory, '-o');
xlabel('Epoch');
ylabel('Validation accuracy (%)');
title('Validation accuracy');
grid on;

%% ---- Save model ---------------------------------------------------------

outModel = fullfile(fileparts(datasetFile), ...
    ['Step5_lowmem_2DCNN_' datestr(now,'yyyymmdd_HHMMSS') '.mat']);

save(outModel, ...
    'net', 'datasetFile', 'trainLossHistory', 'valAccHistory', ...
    'valAcc','valLoss','testAcc','testLoss','yPredVal','yPredTest', ...
    'miniBatchSize','numEpochs','learnRate', '-v7.3');

fprintf('\nSaved model:\n%s\n', outModel);

fprintf('\nStep 5 complete.\n');

%% ========================================================================
% Helper functions
% ========================================================================

function [dlX, dlT] = readBatchFromMatfile(m, xVar, y, idx, C, T, useGPU)

    % Read only this batch from the MAT file.
    switch xVar
        case 'X_train'
            X = m.X_train(idx,:,:);
        case 'X_val'
            X = m.X_val(idx,:,:);
        case 'X_test'
            X = m.X_test(idx,:,:);
        otherwise
            error('Unknown xVar: %s', xVar);
    end

    % X is batch x channels x time.
    B = size(X,1);

    X = reshape(X, [B C T 1]);      % B x C x T x 1
    X = permute(X, [2 3 4 1]);      % C x T x 1 x B

    dlX = dlarray(single(X), 'SSCB');

    yb = y(idx);

    Tmat = zeros(2, B, 'single');
    Tmat(1, yb == 0) = 1;   % NoSpike
    Tmat(2, yb == 1) = 1;   % Spike

    dlT = dlarray(Tmat, 'CB');

    if useGPU
        dlX = gpuArray(dlX);
        dlT = gpuArray(dlT);
    end
end

function [loss, gradients, state] = modelGradients(net, dlX, dlT)

    [dlY, state] = forward(net, dlX);

    % Manual cross-entropy to avoid version-specific options.
    epsVal = 1e-7;
    loss = -mean(sum(dlT .* log(dlY + epsVal), 1));

    gradients = dlgradient(loss, net.Learnables);
end

function [acc, meanLoss, yPredAll] = evaluateModelLowMem(net, m, xVar, y, miniBatchSize, C, T, useGPU)

    N = numel(y);

    yPredAll = zeros(N,1);
    lossSum = 0;
    nBatches = 0;

    for bStart = 1:miniBatchSize:N

        bEnd = min(bStart + miniBatchSize - 1, N);
        idx = bStart:bEnd;

        [dlX, dlT] = readBatchFromMatfile(m, xVar, y, idx, C, T, useGPU);

        dlY = predict(net, dlX);

        epsVal = 1e-7;
        loss = -mean(sum(dlT .* log(dlY + epsVal), 1));

        P = gather(extractdata(dlY));  % 2 x B

        [~, predClass] = max(P, [], 1);
        yPred = predClass(:) - 1;      % convert 1/2 to 0/1

        yPredAll(idx) = yPred;

        lossSum = lossSum + double(gather(extractdata(loss)));
        nBatches = nBatches + 1;
    end

    acc = mean(yPredAll == y);
    meanLoss = lossSum / max(nBatches,1);
end