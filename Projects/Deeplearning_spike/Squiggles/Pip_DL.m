%% The Spike Detection MEG pipline

% Spike Detection MEG pipline, deep learning process
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 07/02/2023

clear; clc, close('all'); warning off

%% FieldTrip toolbox
restoredefaultpath % reset the default path
ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
addpath(ft_path);
ft_defaults

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new/')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/helper')

%%
datadir = '/data/MEG/Research/SpikeDectection/Epil_annotated_data/annotated_data';

cd(datadir)
d_spike = rdir([datadir,'/*.mat']);

%
clear subj run sub_run
for i=1:length(d_spike)
    [pathstr, name] = fileparts(d_spike(i).name);
    tkz = tokenize(name,'_');
    subj{i} = [tkz{1}, '_', tkz{2}];
    sub_run{i,:} = [subj{i}, '_', tkz{3}];
end
[sub_run_unq,IA,IC] = unique(sub_run);
% disp(sub_run_unq);

sub_run_unq_spike = [];
for i=1:length(sub_run_unq)
    sub_run_unq_spike{i} = [num2str(i), '_', sub_run_unq{i}];
end
disp(sub_run_unq_spike')

%%
datadir = '/data/MEG/Research/SpikeDectection/Epil_annotated_data/annotated_data_nospike';
cd(datadir)
d_nospike = rdir([datadir,'/*.mat']);

clear subj run sub_run
for i=1:length(d_nospike)
    [pathstr, name] = fileparts(d_nospike(i).name);
    tkz = tokenize(name,'_');
    subj{i} = [tkz{1}, '_', tkz{2}];
    sub_run{i,:} = [subj{i}, '_', tkz{3}];
end
[sub_run_unq,IA,IC] = unique(sub_run);
% disp(sub_run_unq);

sub_run_unq_nospike = [];
for i=1:length(sub_run_unq)
    sub_run_unq_nospike{i} = [num2str(i), '_', sub_run_unq{i}];
end
disp(sub_run_unq_nospike')

%% BAD DATA?
% baddata = [12];
% ndata = [1:11,13:89,91:length(d_spike)];

%% Spike data
ft_progress('init', 'text',     'please wait ...');

anot_spike_apd = [];
for i= ndata
    ft_progress(i/length(ndata), 'Reading spike data %d from %d', i, length(ndata));
    [pathstr, name] = fileparts(d_spike(i).name);
    spikedata = load(d_spike(i).name);
    anot_spike = [];
    for j=1:length(spikedata.anot_data_all)
        anot_data = spikedata.anot_data_all{j};
        anot_data = spikedata.anot_data_all{j};
        anot_spike(:,:,j) = do_normalize_data(anot_data.trial{1}(:,1:67));
    end
    anot_spike_apd = cat(3, anot_spike_apd, anot_spike);
end
ft_progress('close');

size(anot_spike_apd)

%% non-spike data
ft_progress('init', 'text',     'please wait ...');

anot_nospike_apd = [];
for i= ndata
    ft_progress(i/length(ndata), 'Reading non-spike data %d from %d', i, length(ndata));
    [pathstr, name] = fileparts(d_nospike(i).name);
    no_spikedata = load(d_nospike(i).name);
    anot_nonspike = [];  
    for j=1:length(no_spikedata.anot_data_all)
        anot_data = no_spikedata.anot_data_all{j};
        anot_nonspike(:,:,j) = do_normalize_data (anot_data.trial{1}(:,1:67));
    end
    anot_nospike_apd = cat(3, anot_nospike_apd, anot_nonspike);
end
ft_progress('close');

%%
% Convert data into suitable format for LSTM: cell array of 2D matrices
size(anot_spike_apd)
numTrials = size(anot_spike_apd, 3);
dataCell_spike = cell(numTrials, 1);
for i = 1:numTrials
    dataCell_spike{i} = squeeze(anot_spike_apd(:,:,i));
end
labels_spike = ones(1, length(dataCell_spike));
size(dataCell_spike)

size(anot_nospike_apd)
numTrials = size(anot_nospike_apd, 3);
dataCell_nospike = cell(numTrials, 1);
for i = 1:numTrials
    dataCell_nospike{i} = squeeze(anot_nospike_apd(:,:,i));
end
labels_nospike = zeros(1, length(dataCell_nospike));
size(dataCell_nospike)

%% LSTM net
inputSize = size(anot_spike_apd, 1); % number of channels
numHiddenUnits1 = 100;
numHiddenUnits2 = 50;
numClasses = 2; % spike or not spike

layers = [ ...
    sequenceInputLayer(inputSize)
    lstmLayer(numHiddenUnits1,'OutputMode','sequence')
    dropoutLayer(0.2)
    lstmLayer(numHiddenUnits2,'OutputMode','last')
    dropoutLayer(0.2)
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];

% % Set options for training
options = trainingOptions('adam', ...
    'MaxEpochs',100, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.005, ...
    'MiniBatchSize', 32, ... % Try different batch sizes
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',20, ...
    'LearnRateDropFactor',0.2, ...
    'Verbose',1);

%% removing nan values (bad data!)
containsNaN1 = cellfun(@(x) any(isnan(x(:))), dataCell_spike);
containsNaN2 = cellfun(@(x) any(isnan(x(:))), dataCell_nospike);

idx1  = find(containsNaN1); idx2  = find(containsNaN2);
idx = union(idx1, idx2);

dataCell_spike_new = dataCell_spike; dataCell_spike_new(idx) = [];
dataCell_nospike_new = dataCell_nospike; dataCell_nospike_new(idx) = [];
labels_spike_new = labels_spike; labels_spike_new(idx) = [];
labels_nospike_new = labels_nospike; labels_nospike_new(idx) = [];

%% Train full model
dataCell = [dataCell_spike_new;dataCell_nospike_new];
numericLabels = [labels_spike_new,labels_nospike_new]';

% Convert numeric labels to categorical
categories = {'NoSpike', 'Spike'};
categoricalLabels = categorical(numericLabels, [0 1], categories);

size(categoricalLabels); size(dataCell)

%
% Data normalization (using z-score normalization as an example)
% Assume that dataCell is a cell array of 2D matrices (channels x time)
numTrials = numel(dataCell);
for i = 1:numTrials
    dataCell{i} = (dataCell{i} - mean(dataCell{i}, 'all')) ./ std(dataCell{i}, 0, 'all');
end

% Now 'layers' should be a valid network architecture for training
netCombined = trainNetwork(dataCell, categoricalLabels, layers, options);
% savelebel = size(dataCell,1);
% save(['/data/MEG/Research/awang/Scripts/net_', num2str(savelebel),'.mat'], 'netCombined')

disp('1:yes')
disp('0:no')
snet = input('save network:');

if snet ==1
    savelebel = size(dataCell,1);
    saved_folder = '/data/MEG/Research/SpikeDectection/Scripts/trained model';
    save(fullfile(saved_folder,['/net_', num2str(savelebel),'.mat']), 'netCombined')
end

%%
numSamples  = size(dataCell,1);
randIndices = randperm(numSamples);

subsetIndices = randIndices(1:round(numSamples/3));
XSubset = dataCell(subsetIndices, :);
YSubset = categoricalLabels(subsetIndices);

YPred = classify(netCombined, XSubset);

% Calculate accuracy on the subset
accuracy = sum(YPred == YSubset) / numel(YSubset);
disp(['Accuracy: ', num2str(accuracy)])

%% TEST
numSamples  = size(dataCell,1);
randIndices = randperm(numSamples);

n = 122;

XSubset = dataCell(n, :);
YSubset = categoricalLabels(n);

YPred = classify(netCombined, XSubset);

% Calculate accuracy on the subset
accuracy = sum(YPred == YSubset) / numel(YSubset);
disp(['Accuracy: ', num2str(accuracy)])

%% Cross-validation
labels = categoricalLabels;
k = 10;  % number of folds
indices = crossvalind('Kfold', numel(labels), k);  % create indices for the k folds

accuracies = zeros(1, k);  % vector to store the accuracy for each fold

for i = 1:k
    % Split the data into training and validation sets
    trainData = dataCell(indices ~= i);  % training set is all data not in fold i
    trainLabels = labels(indices ~= i);  % corresponding training labels
    validationData = dataCell(indices == i);  % validation set is data in fold i
    validationLabels = labels(indices == i);  % corresponding validation labels
    
    %     % Define the layers
    %     layers = [ ...
    %         sequenceInputLayer(inputSize)
    %         lstmLayer(numHiddenUnits, 'OutputMode', 'last')
    %         dropoutLayer(0.5)
    %         fullyConnectedLayer(numClasses)
    %         softmaxLayer
    %         classificationLayer];
    %
    %
    %     % Define the options
    %     options = trainingOptions('adam', ...
    %         'MaxEpochs',100, ...
    %         'GradientThreshold',1, ...
    %         'InitialLearnRate',0.005, ...
    %         'LearnRateSchedule','piecewise', ...
    %         'LearnRateDropPeriod',20, ...
    %         'LearnRateDropFactor',0.2, ...
    %         'L2Regularization', 0.001, ...
    %         'Verbose',1);
    
    
    layers = [ ...
        sequenceInputLayer(inputSize)
        lstmLayer(numHiddenUnits1,'OutputMode','sequence')
        dropoutLayer(0.2)
        lstmLayer(numHiddenUnits2,'OutputMode','last')
        dropoutLayer(0.2)
        fullyConnectedLayer(numClasses)
        softmaxLayer
        classificationLayer];
    
    % % Set options for training
    options = trainingOptions('adam', ...
        'MaxEpochs',100, ...
        'GradientThreshold',1, ...
        'InitialLearnRate',0.005, ...
        'MiniBatchSize', 32, ... % Try different batch sizes
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropPeriod',20, ...
        'LearnRateDropFactor',0.2, ...
        'Verbose',1);
    
    % Train the model on the training data
    net = trainNetwork(trainData, categorical(trainLabels), layers, options);
    
    % Evaluate the model on the validation data
    YPred = classify(net, validationData);
    accuracies(i) = sum(YPred == categorical(validationLabels)) / numel(validationLabels);
    disp(accuracies(i))
end

% Average accuracy across all folds
avg_accuracy = mean(accuracies);

%%
% outsum_all = [];
% % for i= subid%1:length(d)
% for i= 1:20%length(d)
%     disp([num2str(i),'/',num2str(length(d))])
%     [pathstr, name] = fileparts(d(i).name);
%     load(d(i).name);
%     
%     outsum = [];
%     for j=1:length(anot_data_all)
%         anot_data = anot_data_all{j};
%         
%         cfg = [];
%         cfg.plot = 0;
%         outsum(j,:,:) = do_conn(cfg,anot_data.trial{1});
%         
%         cfg = [];
%         cfg.blocksize = anot_data.time{1}(end) - anot_data.time{1}(1);
%         cfg.viewmode = 'vertical'; %butterfly';
%         cfg.continuous = 'yes';
%         cfg.axisfontsize = 7;
%         cfg.fontsize = 7;
%         cfg.channel = 'EEG*';
%         cfg.preproc.demean = 'yes';
%         cfg.position = [300   900   500   1500];
%         ft_databrowser(cfg, anot_data);
%         cfg.channel = 'MEG*';
%         cfg.position = [850   900   500   1500];
%         ft_databrowser(cfg, anot_data);
%         
%         cfg = [];
%         cfg.channel = 'MEG*';
%         MEG_data = ft_selectdata(cfg, anot_data);
%         cfg.channel = 'EEG*';
%         EEG_data = ft_selectdata(cfg, anot_data);
%         
%         tmp = rms(anot_data.trial{1});
%         tmp = rms(MEG_data.trial{1});
%         tmp = kurtosis(MEG_data.trial{1});
%         tmp = kurtosis(MEG_data.trial{1}) ./rms(MEG_data.trial{1});
%         tmp = kurtosis(EEG_data.trial{1}) ./rms(EEG_data.trial{1});
%         tmp = kurtosis(EEG_data.trial{1});
%         figure, plot(anot_data.time{1},tmp)
%         %
%         %         pause,
%         %         close all,
%     end
%     outsum_all{i} = squeeze(mean(outsum,1));
% end
