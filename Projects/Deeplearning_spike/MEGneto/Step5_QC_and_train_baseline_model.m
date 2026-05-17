%% ------------------------------------------------------------------------
% Step5_QC_and_train_baseline_model.m
%
% Load dataset_DL_ready.mat, check class balance and dimensions,
% then train a simple baseline classifier.
% -------------------------------------------------------------------------

clc; clear;

datasetFile = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files/dataset_DL_ready.mat';

%% ---- Inspect saved variables without loading all data -------------------

whos('-file', datasetFile)

m = matfile(datasetFile);

fprintf('\nArray sizes:\n');
disp(size(m, 'X_train'))
disp(size(m, 'X_val'))
disp(size(m, 'X_test'))

S = load(datasetFile, ...
    'dataset', ...
    'y_train','y_val','y_test', ...
    'subject_id_train','subject_id_val','subject_id_test', ...
    'base_train','base_val','base_test', ...
    'Tsummary');

dataset = S.dataset;

y_train = S.y_train;
y_val   = S.y_val;
y_test  = S.y_test;

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

%% ---- Load arrays for first baseline model -------------------------------
% If this causes memory issue, we can switch this part to minibatch/datastore.

X_train = m.X_train;
X_val   = m.X_val;
X_test  = m.X_test;

% X is N x C x T = N x 306 x 306
% MATLAB sequenceInputLayer expects cell array, each cell = features x time.

XTrainSeq = squeeze_to_sequence_cells(X_train);
XValSeq   = squeeze_to_sequence_cells(X_val);
XTestSeq  = squeeze_to_sequence_cells(X_test);

YTrain = categorical(y_train(:), [0 1], {'NoSpike','Spike'});
YVal   = categorical(y_val(:),   [0 1], {'NoSpike','Spike'});
YTest  = categorical(y_test(:),  [0 1], {'NoSpike','Spike'});

%% ---- Simple baseline CNN over time --------------------------------------

numChannels = size(X_train,2);

layers = [
    sequenceInputLayer(numChannels, 'Name','input')

    convolution1dLayer(7, 32, 'Padding','same', 'Name','conv1')
    batchNormalizationLayer('Name','bn1')
    reluLayer('Name','relu1')
    maxPooling1dLayer(2, 'Stride',2, 'Name','pool1')

    convolution1dLayer(5, 64, 'Padding','same', 'Name','conv2')
    batchNormalizationLayer('Name','bn2')
    reluLayer('Name','relu2')
    maxPooling1dLayer(2, 'Stride',2, 'Name','pool2')

    convolution1dLayer(3, 128, 'Padding','same', 'Name','conv3')
    batchNormalizationLayer('Name','bn3')
    reluLayer('Name','relu3')

    globalAveragePooling1dLayer('Name','gap')

    dropoutLayer(0.3, 'Name','dropout')
    fullyConnectedLayer(2, 'Name','fc')
    softmaxLayer('Name','softmax')
    classificationLayer('Name','classOutput')
];

miniBatchSize = 64;

options = trainingOptions('adam', ...
    'MaxEpochs', 20, ...
    'MiniBatchSize', miniBatchSize, ...
    'InitialLearnRate', 1e-3, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{XValSeq,YVal}, ...
    'ValidationFrequency', 50, ...
    'Verbose', true, ...
    'Plots','training-progress');

net = trainNetwork(XTrainSeq, YTrain, layers, options);

%% ---- Evaluate -----------------------------------------------------------

YPredVal  = classify(net, XValSeq,  'MiniBatchSize', miniBatchSize);
YPredTest = classify(net, XTestSeq, 'MiniBatchSize', miniBatchSize);

accVal  = mean(YPredVal  == YVal);
accTest = mean(YPredTest == YTest);

fprintf('\nValidation accuracy: %.2f%%\n', 100*accVal);
fprintf('Test accuracy:       %.2f%%\n', 100*accTest);

figure;
confusionchart(YTest, YPredTest);
title('Test confusion matrix');

%% ---- Save model ---------------------------------------------------------

outModel = fullfile(fileparts(datasetFile), ...
    ['Step5_baseline_CNN1D_' datestr(now,'yyyymmdd_HHMMSS') '.mat']);

save(outModel, 'net', 'options', 'accVal', 'accTest', 'datasetFile', '-v7.3');

fprintf('\nSaved model:\n%s\n', outModel);

%% ======================== helper function ================================

function Xseq = squeeze_to_sequence_cells(X)

    % Input X: N x C x T
    N = size(X,1);

    Xseq = cell(N,1);

    for i = 1:N
        Xseq{i} = squeeze(X(i,:,:));  % C x T
    end
end