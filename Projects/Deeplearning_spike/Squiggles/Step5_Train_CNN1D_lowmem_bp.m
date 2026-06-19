%% ------------------------------------------------------------------------
% Step5_Train_CNN1D_lowmem_v3.m
%
% Low-memory CNN1D training from dataset_DL_ready*.mat.
%
% Supports:
%   Format A: top-level X_train/X_val/X_test variables
%   Format B: dataset.X_train/X_val/X_test inside dataset struct
%
% Input data format:
%   X: N x channels x time
%
% Network input format:
%   dlarray C x T x B, label 'CTB'
%
% Main improvements:
%   - CNN1D architecture instead of stacked LSTM
%   - contiguous balanced batch reading
%   - best model saved as bestNet for Step 6 compatibility
%   - early stopping using validation loss
%   - gradient clipping
%   - class-weighted loss, no-op when classes are balanced
% -------------------------------------------------------------------------

clc; clear;

%% ======================== USER SETTINGS =================================

% Full-window dataset:
% datasetFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Squiggles/dataset_DL_ready.mat';

% Narrow-window dataset:
% datasetFile = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Squiggles/dataset_DL_ready_win_m250_p500.mat';
datasetFile = '/home/vyoussofzadeh/Data/DL_model/dataset_DL_ready_win_m250_p500_bp5_50.mat';
% datasetFile = '/home/vyoussofzadeh/Data/DL_model/dataset_DL_ready_win_m250_p500.mat';

modelTag = 'CNN1D_bp5_50_m250_p500';

% Label convention
LABEL_NOSPIKE = 0;
LABEL_SPIKE   = 1;

% Training settings
numEpochs     = 15;
halfBatch     = 32;          % 32 spike + 32 no-spike = 64 total
miniBatchSize = 2 * halfBatch;
learnRate     = 2e-4;
gradClipVal   = 1.0;

% Early stopping based on validation loss
patienceLimit = 4;
minDelta      = 1e-4;

% Evaluation batch size
evalMiniBatchSize = 128;

useGPU = canUseGPU();
fprintf('\nUse GPU: %d\n', useGPU);

rng(1);

%% ======================== OPEN DATASET ==================================

fprintf('\nDataset file:\n%s\n', datasetFile);
whos('-file', datasetFile)

[data, dataset, y_train, y_val, y_test, subject_id_train, subject_id_val, subject_id_test] = ...
    openDataset(datasetFile);

szTrain = getXSize(data, 'X_train');
szVal   = getXSize(data, 'X_val');
szTest  = getXSize(data, 'X_test');

fprintf('\nArray sizes:\n');
fprintf('  X_train : [%s]\n', num2str(szTrain));
fprintf('  X_val   : [%s]\n', num2str(szVal));
fprintf('  X_test  : [%s]\n', num2str(szTest));

fprintf('\nClass counts:\n');
fprintf('  Train : spike=%d, nospike=%d\n', sum(y_train==LABEL_SPIKE), sum(y_train==LABEL_NOSPIKE));
fprintf('  Val   : spike=%d, nospike=%d\n', sum(y_val  ==LABEL_SPIKE), sum(y_val  ==LABEL_NOSPIKE));
fprintf('  Test  : spike=%d, nospike=%d\n', sum(y_test ==LABEL_SPIKE), sum(y_test ==LABEL_NOSPIKE));

fprintf('\nSubject leakage check:\n');
trainSubj = unique(subject_id_train);
valSubj   = unique(subject_id_val);
testSubj  = unique(subject_id_test);

fprintf('  Train-Val overlap  : %d\n', numel(intersect(trainSubj, valSubj)));
fprintf('  Train-Test overlap : %d\n', numel(intersect(trainSubj, testSubj)));
fprintf('  Val-Test overlap   : %d\n', numel(intersect(valSubj, testSubj)));

Ntrain = szTrain(1);
C      = szTrain(2);
T      = szTrain(3);

numClasses = 2;

%% ======================== CLASS WEIGHTS =================================

nPos = sum(y_train == LABEL_SPIKE);
nNeg = sum(y_train == LABEL_NOSPIKE);

wPos = (nPos + nNeg) / (2 * nPos);
wNeg = (nPos + nNeg) / (2 * nNeg);

classWeights = single([wNeg; wPos]);   % row 1 = NoSpike, row 2 = Spike

fprintf('\nClass weights: NoSpike=%.3f  Spike=%.3f\n', wNeg, wPos);

%% ======================== NETWORK ARCHITECTURE ===========================

layers = [
    sequenceInputLayer(C, 'Normalization','none', 'Name','input', 'MinLength', T)

    convolution1dLayer(7, 64, 'Padding','same', 'Name','conv1')
    batchNormalizationLayer('Name','bn1')
    reluLayer('Name','relu1')
    maxPooling1dLayer(2, 'Stride',2, 'Name','pool1')

    convolution1dLayer(5, 128, 'Padding','same', 'Name','conv2')
    batchNormalizationLayer('Name','bn2')
    reluLayer('Name','relu2')
    maxPooling1dLayer(2, 'Stride',2, 'Name','pool2')

    convolution1dLayer(3, 128, 'Padding','same', 'Name','conv3')
    batchNormalizationLayer('Name','bn3')
    reluLayer('Name','relu3')

    globalAveragePooling1dLayer('Name','gap')

    dropoutLayer(0.4, 'Name','dropout')
    fullyConnectedLayer(numClasses, 'Name','fc')
    softmaxLayer('Name','softmax')
];

net = dlnetwork(layers);

disp(net)

%% ======================== BALANCED TRAIN BLOCKS ==========================

posIdx = find(y_train == LABEL_SPIKE);
negIdx = find(y_train == LABEL_NOSPIKE);

fprintf('\nTrain rows: pos=%d, neg=%d\n', numel(posIdx), numel(negIdx));
fprintf('  First/last pos rows: %d - %d\n', posIdx(1), posIdx(end));
fprintf('  First/last neg rows: %d - %d\n', negIdx(1), negIdx(end));

posBlocks = makeContigBlocks(posIdx, halfBatch);
negBlocks = makeContigBlocks(negIdx, halfBatch);

% Drop partial blocks so every training batch is same size
posBlocks = posBlocks(cellfun(@numel, posBlocks) == halfBatch);
negBlocks = negBlocks(cellfun(@numel, negBlocks) == halfBatch);

nTrainBlocks = min(numel(posBlocks), numel(negBlocks));

if nTrainBlocks == 0
    error('No balanced training blocks created. Check halfBatch and class counts.');
end

fprintf('Balanced train blocks per epoch: %d\n', nTrainBlocks);

printEvery = max(1, min(25, floor(nTrainBlocks / 10)));

%% ======================== INITIALIZE HISTORY =============================

trainLossHistory = zeros(numEpochs, 1);
valAccHistory    = zeros(numEpochs, 1);
valLossHistory   = zeros(numEpochs, 1);

bestNet     = net;
netFinal    = net;
bestEpoch   = 0;
bestValAcc  = -inf;
bestValLoss = inf;
patience    = 0;

trailingAvg   = [];
trailingAvgSq = [];
iteration     = 0;

%% ======================== TRAINING LOOP ==================================

fprintf('\nStarting CNN1D low-memory training...\n');
fprintf('numEpochs=%d | halfBatch=%d | miniBatchSize=%d | learnRate=%g\n\n', ...
    numEpochs, halfBatch, miniBatchSize, learnRate);

actualEpochs = numEpochs;

for epoch = 1:numEpochs

    tEpoch = tic;

    posOrder = randperm(numel(posBlocks));
    negOrder = randperm(numel(negBlocks));

    epochLoss = 0;

    for bb = 1:nTrainBlocks

        iteration = iteration + 1;

        posRows = posBlocks{posOrder(bb)};
        negRows = negBlocks{negOrder(bb)};

        [dlX, dlT] = readBalancedBatch(data, y_train, posRows, negRows, C, T, useGPU, LABEL_SPIKE);

        [loss, gradients, state] = dlfeval(@modelGradients, net, dlX, dlT, classWeights);

        % Update batch norm state
        net.State = state;

        % Per-parameter gradient clipping
        gradients = dlupdate(@(g) clipGradient(g, gradClipVal), gradients);

        [net, trailingAvg, trailingAvgSq] = adamupdate( ...
            net, gradients, trailingAvg, trailingAvgSq, iteration, learnRate);

        epochLoss = epochLoss + double(gather(extractdata(loss)));

        if mod(bb, printEvery) == 0
            fprintf('  Epoch %d/%d | batch %d/%d | loss %.4f\n', ...
                epoch, numEpochs, bb, nTrainBlocks, ...
                double(gather(extractdata(loss))));
        end
    end

    meanTrainLoss = epochLoss / nTrainBlocks;

    fprintf('Evaluating validation...\n');

    [valAcc, valLoss] = evaluateModel( ...
        net, data, 'X_val', y_val, evalMiniBatchSize, C, T, useGPU, LABEL_SPIKE);

    trainLossHistory(epoch) = meanTrainLoss;
    valAccHistory(epoch)    = valAcc;
    valLossHistory(epoch)   = valLoss;

    elapsed = toc(tEpoch);
    eta = elapsed * (numEpochs - epoch);

    fprintf('\nEpoch %d/%d | train loss %.4f | val loss %.4f | val acc %.2f%% | %.1fs elapsed | ETA %.0fs\n', ...
        epoch, numEpochs, meanTrainLoss, valLoss, 100*valAcc, elapsed, eta);

    % Early stopping by validation loss
    if valLoss < bestValLoss - minDelta

        bestValAcc  = valAcc;
        bestValLoss = valLoss;
        bestNet     = net;
        bestEpoch   = epoch;
        patience    = 0;

        fprintf('  *** New best model: val acc %.2f%%, val loss %.4f ***\n', ...
            100*valAcc, valLoss);

    else

        patience = patience + 1;
        fprintf('  No validation-loss improvement (%d/%d)\n', patience, patienceLimit);
    end

    fprintf('\n');

    if patience >= patienceLimit
        fprintf('Early stopping triggered at epoch %d.\n', epoch);
        actualEpochs = epoch;
        break
    end
end

netFinal = net;

trainLossHistory = trainLossHistory(1:actualEpochs);
valAccHistory    = valAccHistory(1:actualEpochs);
valLossHistory   = valLossHistory(1:actualEpochs);

%% ======================== FINAL EVALUATION ===============================

net = bestNet;

fprintf('\nFinal evaluation using best model from epoch %d...\n', bestEpoch);

[valAcc, valLoss, yPredVal, pVal] = evaluateModel( ...
    net, data, 'X_val', y_val, evalMiniBatchSize, C, T, useGPU, LABEL_SPIKE);

[testAcc, testLoss, yPredTest, pTest] = evaluateModel( ...
    net, data, 'X_test', y_test, evalMiniBatchSize, C, T, useGPU, LABEL_SPIKE);

Mval  = metrics_binary(y_val,  yPredVal,  pVal);
Mtest = metrics_binary(y_test, yPredTest, pTest);

fprintf('\n--- Final results: %s, best epoch %d ---\n', modelTag, bestEpoch);
print_metrics('Val', Mval);
print_metrics('Test', Mtest);

fprintf('\nTest confusion matrix:\n');
fprintf('                 Pred NoSpike    Pred Spike\n');
fprintf('True NoSpike      %8d      %8d\n', Mtest.TN, Mtest.FP);
fprintf('True Spike        %8d      %8d\n', Mtest.FN, Mtest.TP);

%% ======================== PLOTS ==========================================

figure('Name','Test confusion matrix');
confusionchart( ...
    categorical(y_test,    [LABEL_NOSPIKE LABEL_SPIKE], {'NoSpike','Spike'}), ...
    categorical(yPredTest, [LABEL_NOSPIKE LABEL_SPIKE], {'NoSpike','Spike'}));
title(sprintf('Test confusion matrix: %s', modelTag));

figure('Name','Training curves');
tiledlayout(2,1);

nexttile;
plot(trainLossHistory, '-o', 'DisplayName','Train loss');
hold on;
plot(valLossHistory, '-s', 'DisplayName','Val loss');
xlabel('Epoch');
ylabel('Loss');
title('Loss');
legend;
grid on;

nexttile;
plot(100*valAccHistory, '-o');
xlabel('Epoch');
ylabel('Validation accuracy (%)');
title('Validation accuracy');
grid on;

%% ======================== SAVE MODEL =====================================

timestamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));

outModel = fullfile(fileparts(datasetFile), ...
    sprintf('Step5_lowmem_%s_%s.mat', modelTag, timestamp));

save(outModel, ...
    'bestNet','netFinal','bestEpoch', ...
    'datasetFile','modelTag', ...
    'trainLossHistory','valAccHistory','valLossHistory', ...
    'valAcc','valLoss','testAcc','testLoss', ...
    'yPredVal','yPredTest','pVal','pTest', ...
    'Mval','Mtest', ...
    'miniBatchSize','halfBatch','numEpochs','actualEpochs','learnRate','gradClipVal', ...
    'classWeights','LABEL_NOSPIKE','LABEL_SPIKE', ...
    '-v7.3');

fprintf('\nSaved model:\n%s\n', outModel);
fprintf('\nStep 5 complete.\n');

%% ========================================================================
% Helper functions
% ========================================================================

function [data, dataset, y_train, y_val, y_test, subject_id_train, subject_id_val, subject_id_test] = openDataset(datasetFile)

    info = whos('-file', datasetFile);
    names = {info.name};

    if all(ismember({'X_train','X_val','X_test','y_train','y_val','y_test'}, names))

        data.mode = 'matfile';
        data.m = matfile(datasetFile);

        S = load(datasetFile, ...
            'dataset', ...
            'y_train','y_val','y_test', ...
            'subject_id_train','subject_id_val','subject_id_test');

        if isfield(S, 'dataset')
            dataset = S.dataset;
        else
            dataset = struct();
        end

        y_train = double(S.y_train(:));
        y_val   = double(S.y_val(:));
        y_test  = double(S.y_test(:));

        subject_id_train = S.subject_id_train;
        subject_id_val   = S.subject_id_val;
        subject_id_test  = S.subject_id_test;

    elseif ismember('dataset', names)

        S = load(datasetFile, 'dataset');
        dataset = S.dataset;

        if ~isfield(dataset,'X_train') || ~isfield(dataset,'X_val') || ~isfield(dataset,'X_test')
            error('dataset exists but does not contain X_train/X_val/X_test.');
        end

        data.mode = 'struct';
        data.dataset = dataset;

        y_train = double(dataset.y_train(:));
        y_val   = double(dataset.y_val(:));
        y_test  = double(dataset.y_test(:));

        subject_id_train = dataset.subject_id_train;
        subject_id_val   = dataset.subject_id_val;
        subject_id_test  = dataset.subject_id_test;

    else
        error('Unknown dataset format: %s', datasetFile);
    end
end

function sz = getXSize(data, xVar)

    switch data.mode

        case 'matfile'
            sz = size(data.m, xVar);

        case 'struct'
            sz = size(data.dataset.(xVar));

        otherwise
            error('Unknown data mode: %s', data.mode);
    end
end

function [dlX, dlT] = readBalancedBatch(data, y_train, posRows, negRows, C, T, useGPU, LABEL_SPIKE)

    Xpos = readX(data, 'X_train', posRows, C, T);
    Xneg = readX(data, 'X_train', negRows, C, T);

    X = single(cat(1, Xpos, Xneg));
    yb = y_train([posRows(:); negRows(:)]);

    B = size(X,1);

    ord = randperm(B);
    X = X(ord,:,:);
    yb = yb(ord);

    dlX = dlarray(permute(X, [2 3 1]), 'CTB');

    Tmat = zeros(2, B, 'single');
    Tmat(1, yb ~= LABEL_SPIKE) = 1;
    Tmat(2, yb == LABEL_SPIKE) = 1;

    dlT = dlarray(Tmat, 'CB');

    if useGPU
        dlX = gpuArray(dlX);
        dlT = gpuArray(dlT);
    end
end

function X = readX(data, xVar, idx, C, T)

    idx = idx(:);

    switch data.mode

        case 'matfile'

            % If contiguous, read one block. Otherwise read safely.
            if numel(idx) == idx(end) - idx(1) + 1
                switch xVar
                    case 'X_train'
                        X = data.m.X_train(idx(1):idx(end),:,:);
                    case 'X_val'
                        X = data.m.X_val(idx(1):idx(end),:,:);
                    case 'X_test'
                        X = data.m.X_test(idx(1):idx(end),:,:);
                    otherwise
                        error('Unknown xVar: %s', xVar);
                end
            else
                X = readNonContiguousMatfile(data.m, xVar, idx, C, T);
            end

        case 'struct'

            X = data.dataset.(xVar)(idx,:,:);

        otherwise
            error('Unknown data mode: %s', data.mode);
    end

    if size(X,1) ~= numel(idx)
        X = reshape(X, [numel(idx) C T]);
    end

    X = single(X);
end

function X = readNonContiguousMatfile(m, xVar, idx, C, T)

    idx = idx(:);
    B = numel(idx);

    [idxSorted, sortOrder] = sort(idx);
    Xsorted = zeros(B, C, T, 'single');

    breaks = [1; find(diff(idxSorted) > 1) + 1; B + 1];

    writePos = 1;

    for r = 1:numel(breaks)-1

        a = breaks(r);
        b = breaks(r+1) - 1;

        rowStart = idxSorted(a);
        rowEnd   = idxSorted(b);
        nRows    = b - a + 1;

        switch xVar
            case 'X_train'
                Xpart = m.X_train(rowStart:rowEnd,:,:);
            case 'X_val'
                Xpart = m.X_val(rowStart:rowEnd,:,:);
            case 'X_test'
                Xpart = m.X_test(rowStart:rowEnd,:,:);
            otherwise
                error('Unknown xVar: %s', xVar);
        end

        Xsorted(writePos:writePos+nRows-1,:,:) = single(reshape(Xpart, [nRows C T]));
        writePos = writePos + nRows;
    end

    X = zeros(B, C, T, 'single');
    X(sortOrder,:,:) = Xsorted;
end

function [loss, gradients, state] = modelGradients(net, dlX, dlT, classWeights)

    [dlY, state] = forward(net, dlX);

    w = classWeights(1) .* dlT(1,:) + classWeights(2) .* dlT(2,:);

    epsVal = 1e-7;
    loss = -mean(w .* sum(dlT .* log(dlY + epsVal), 1));

    gradients = dlgradient(loss, net.Learnables);
end

function g = clipGradient(g, threshold)

    if isempty(g)
        return;
    end

    gData = extractdata(g);
    nrm = sqrt(sum(gData.^2, 'all'));

    if nrm > threshold
        g = g .* (threshold ./ nrm);
    end
end

function [acc, meanLoss, yPredAll, pSpikeAll] = evaluateModel( ...
    net, data, xVar, y, miniBatchSize, C, T, useGPU, LABEL_SPIKE)

    N = numel(y);

    yPredAll = zeros(N,1);
    pSpikeAll = zeros(N,1);

    lossSum = 0;
    nBatches = 0;

    for bStart = 1:miniBatchSize:N

        bEnd = min(bStart + miniBatchSize - 1, N);
        idx = (bStart:bEnd)';

        X = readX(data, xVar, idx, C, T);

        [dlX, dlT] = makeDlBatch(net, X, y(idx), useGPU, LABEL_SPIKE);

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

function [dlX, dlT] = makeDlBatch(net, X, y, useGPU, LABEL_SPIKE)

    inputLayer = net.Layers(1);
    inputClass = class(inputLayer);

    B = size(X,1);

    if contains(inputClass, 'ImageInputLayer')

        X = reshape(single(X), [B size(X,2) size(X,3) 1]);
        X = permute(X, [2 3 4 1]);
        dlX = dlarray(X, 'SSCB');

    elseif contains(inputClass, 'SequenceInputLayer')

        X = permute(single(X), [2 3 1]);
        dlX = dlarray(X, 'CTB');

    else
        error('Unsupported input layer type: %s', inputClass);
    end

    Tmat = zeros(2, B, 'single');
    Tmat(1, y ~= LABEL_SPIKE) = 1;
    Tmat(2, y == LABEL_SPIKE) = 1;

    dlT = dlarray(Tmat, 'CB');

    if useGPU
        dlX = gpuArray(dlX);
        dlT = gpuArray(dlT);
    end
end

function blocks = makeContigBlocks(idx, blockSize)

    idx = sort(idx(:));
    N = numel(idx);

    blocks = {};

    runStart = 1;
    breakPts = [find(diff(idx) > 1); N];

    for r = 1:numel(breakPts)

        if runStart > N
            break;
        end

        runEnd = breakPts(r);
        thisRun = idx(runStart:runEnd);

        for s = 1:blockSize:numel(thisRun)
            e = min(s + blockSize - 1, numel(thisRun));
            blocks{end+1,1} = thisRun(s:e); %#ok<AGROW>
        end

        runStart = runEnd + 1;
    end
end

function M = metrics_binary(y, pred, prob)

    y = double(y(:));
    pred = double(pred(:));
    prob = double(prob(:));

    TP = sum(y==1 & pred==1);
    TN = sum(y==0 & pred==0);
    FP = sum(y==0 & pred==1);
    FN = sum(y==1 & pred==0);

    acc = mean(pred == y);

    sens = TP / max(TP + FN, 1);
    spec = TN / max(TN + FP, 1);
    balAcc = 0.5 * (sens + spec);

    loss = -mean(y .* log(prob + 1e-7) + (1-y) .* log(1-prob + 1e-7));
    auc = auc_rank(prob, y);

    M = struct();
    M.TP = TP;
    M.TN = TN;
    M.FP = FP;
    M.FN = FN;
    M.acc = acc;
    M.sensitivity = sens;
    M.specificity = spec;
    M.balancedAccuracy = balAcc;
    M.loss = loss;
    M.auc = auc;
end

function auc = auc_rank(scores, labels)

    scores = double(scores(:));
    labels = double(labels(:));

    nPos = sum(labels==1);
    nNeg = sum(labels==0);

    if nPos == 0 || nNeg == 0
        auc = NaN;
        return;
    end

    [scoresSorted, order] = sort(scores);
    ranks = zeros(size(scores));

    i = 1;

    while i <= numel(scoresSorted)

        j = i;

        while j < numel(scoresSorted) && scoresSorted(j+1) == scoresSorted(i)
            j = j + 1;
        end

        avgRank = mean(i:j);
        ranks(order(i:j)) = avgRank;

        i = j + 1;
    end

    sumRanksPos = sum(ranks(labels==1));
    auc = (sumRanksPos - nPos*(nPos+1)/2) / (nPos*nNeg);
end

function print_metrics(name, M)

    fprintf('%s: acc=%.2f%% | balAcc=%.2f%% | sens=%.2f%% | spec=%.2f%% | AUC=%.3f | loss=%.4f\n', ...
        name, 100*M.acc, 100*M.balancedAccuracy, ...
        100*M.sensitivity, 100*M.specificity, M.auc, M.loss);
end