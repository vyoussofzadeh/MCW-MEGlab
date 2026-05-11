function run_spike_classifier_withTESTdata4(varargin)
%% Spike-Detection MEG Pipeline  Spike vs Non-Spike Classifier (Applied Suggestions)
% Deep-learning workflow to classify inter-ictal MEG epochs as **Spike** or **NoSpike**.
% Author: MCW MEG Lab  V. Youssofzadeh <vyoussofzadeh@mcw.edu>
% -------------------------------------------------------------------------
% Usage:
%   run_spike_classifier('arch','transformer','norm','robust','augment',true)
%   run_spike_classifier('arch','lstm','valFrac',0.2)
% -------------------------------------------------------------------------

%% --------------------------- USER SETTINGS ------------------------------
opts.ftRoot      = '/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
opts.codeRoot    = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git';
opts.spikeDir    = '/data/MEG/Research/SpikeDetection/Epil_annotated_data/annotated_data';
opts.noSpikeDir  = '/data/MEG/Research/SpikeDetection/Epil_annotated_data/annotated_data_nospike';
opts.skipList    = {};                 % e.g. {'pilot','bad_subj'}
opts.useGPU      = false;              % auto-downgrades if incompatible
opts.modelDir    = fullfile(opts.spikeDir,'models_classifier');
opts.savePlots   = true;
opts.verbose     = true;
opts.valFrac     = 0.15;
opts.arch        = 'transformer';      % 'lstm' | 'gru_td' | 'transformer'
opts.dModel      = 128;                % transformer width (ignored by LSTM/GRU)
opts.numHeads    = 8;                  % transformer heads (ignored by LSTM/GRU)
opts.norm        = 'zscore';           % 'zscore' | 'robust'
opts.augment     = false;              % light on-the-fly augmentation (train only)
opts.seed        = 0;                  % RNG seed for reproducibility
rng(opts.seed,'twister');

% Allow name-value overrides
opts = parse_inputs(opts,varargin{:});

%% ----------------------------- ENV SETUP --------------------------------
setup_environment(opts);

% GPU detection (explicit)
gpuOK = false;
if opts.useGPU
    try
        g = gpuDevice;
        gpuOK = ~isempty(g) && g.SupportsDouble;
        if gpuOK && opts.verbose
            fprintf('Using GPU: %s (CC %s)\n', g.Name, g.ComputeCapability);
        end
    catch
        gpuOK = false;
    end
end
if ~gpuOK, opts.useGPU = false; disp('Training on CPU (GPU unavailable/incompatible).'); end
if ~exist(opts.modelDir,'dir'); mkdir(opts.modelDir); end

%% ------------------------- LOAD MAT FILE LISTS --------------------------
filesS = list_mat_files(opts.spikeDir);
filesN = list_mat_files(opts.noSpikeDir);
assert(~isempty(filesS) && ~isempty(filesN),'No .mat files found in spike/no-spike directories.');
%
% % Optional skip list
isSkip = @(f) any(cellfun(@(s) contains(f,s), opts.skipList));
filesS = filesS(~arrayfun(@(f) isSkip(f.name), filesS));
filesN = filesN(~arrayfun(@(f) isSkip(f.name), filesN));

%% ----------------------- BUILD TRIAL CELL ARRAYS ------------------------
% trials = build_trials_from_mat(filesS(1), filesN(1));
% trials = build_trials_from_mat(filesS, filesN);
% 
% % Use these in the rest of your pipeline
% data   = trials.data;        % {N×1} cell, each C×T
% labels = trials.labels;      % categorical {'NoSpike','Spike'}
% groups = trials.groups;      % {N×1} group IDs (leakage-safe split key)
% refLbl = trials.refLbl;      % channel labels (row order of data epochs)
% 
% fprintf('Built trials: N=%d  |  C=%d  |  Spike=%d  |  NoSpike=%d\n', ...
%     numel(data), numel(refLbl), sum(labels=='Spike'), sum(labels=='NoSpike'));

load('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/Squiggles/func/DD.mat');

trials = [];
trials.data = data;
trials.labels = labels;
trials.groups = groups;

% fprintf('Built trials: N=%d  |  C=%d  |  Spike=%d  |  NoSpike=%d\n', ...
%     numel(data), numel(refLbl), sum(labels=='Spike'), sum(labels=='NoSpike'));

%% ------------------------ GROUPED TRAIN / VAL SPLIT ---------------------
[trainMask, valMask] = group_holdout_split(groups, opts.valFrac, opts.seed);
XTrain = data(trainMask);   YTrain = labels(trainMask);
XVal   = data(valMask);     YVal   = labels(valMask);

% Keep class order stable across Train/Val
cats   = categories(YTrain);
YTrain = reordercats(YTrain,cats);
YVal   = reordercats(YVal,cats);

assert(~isempty(XTrain), 'Training set is empty after group split. Reduce valFrac or check group IDs.');


%% Smoothing
% After split, before normalization:
Fs = 500;  % <-- set your true sampling rate

% Example A: SavitzkyGolay (very good for preserving spike peaks)
cfg = struct('type','sgolay','frame_ms',9,'order',3);

% Example B: simple symmetric moving average (fast; ok for quick tests)
% cfg = struct('type','movmean','win_ms',7);

% Example C: mild low-pass (overlaps what you already tried)
% cfg = struct('type','lowpass','fc',80,'order',4);

[XTrain, XVal] = smooth_train_val(XTrain, XVal, Fs, cfg);

% ...then proceed to fit_channel_stats(XTrain) and build XTrain2/XVal2

%% ----------------------- CHANNEL-WISE NORMALIZATION ---------------------
switch lower(opts.norm)
    case 'zscore'
        [mu, sd] = fit_channel_stats(XTrain);
        mu = ensure_column(mu); sd = ensure_column(sd);
        zscoref  = @(x) bsxfun(@rdivide, bsxfun(@minus,x,mu), max(sd,1e-12));
        normStats = struct('type','zscore','mu',mu,'sd',sd);
    case 'robust'
        [med, madv] = fit_channel_stats_robust(XTrain);
        med = ensure_column(med); madv = ensure_column(madv);
        zscoref  = @(x) bsxfun(@rdivide, bsxfun(@minus,x,med), max(madv,1e-12));
        normStats = struct('type','robust','med',med,'madv',madv);
        mu = []; sd = [];
    otherwise
        error('Unknown opts.norm: %s', opts.norm);
end
XTrain2 = cellfun(zscoref, XTrain,'UniformOutput',false);
XVal2   = cellfun(zscoref, XVal  ,'UniformOutput',false);

%% ----- Fix a single SequenceLength for both sets (even when PCA is off)
Ttrain = cellfun(@(x) size(x,2), XTrain2);
Tval   = cellfun(@(x) size(x,2), XVal2);
Tfix   = max(1, round(prctile([Ttrain(:); Tval(:)], 95)));
XTrain2 = pad_sequences_to_length(XTrain2, Tfix);
XVal2   = pad_sequences_to_length(XVal2,   Tfix);
padArgs = {"SequenceLength", Tfix, "SequencePaddingDirection","right"};

%% ----- Class weights (use instead of oversampling)
cnt  = countcats(YTrain);
w    = sum(cnt)./max(cnt,1); w = w/mean(w);
clsLayer = classificationLayer('Name','cls','Classes',cats,'ClassWeights',w);

%% ------------------------- BUILD MODEL & OPTIONS ------------------------
% Sanity & input size AFTER any PE addition
chCounts = cellfun(@(x) size(x,1), XTrain2);
assert(~isempty(chCounts),'No training epochs left after cleaning.');
uCh = unique(chCounts);
assert(numel(uCh)==1, 'Inconsistent channel counts across epochs: %s', mat2str(uCh));
inputSize = uCh(1);

% (No second padArgs here  keep the fixed Tfix one)
fprintf('Groups: train=%d | val=%d | Tfix=%d | inputSize=%d\n', ...
    numel(unique(groups(trainMask))), numel(unique(groups(valMask))), Tfix, inputSize);
disp('Class weights:'), disp(array2table(w','VariableNames',cats))

opts.pca = 0;
if opts.pca == 1
    %% ----- DATA REDUCTION (TRAIN ONLY fit, then apply to VAL) -----

    % After normalization + fixed T (XTrain2, YTrain ready)
% fs = 1000;        % <-- set your sampling rate (Hz); or [] for samples
% [Wtrain, compTS] = pca_timeseries_by_class(XTrain2, YTrain, 306, fs);

    
    % 1) Optional: keep only gradiometers (fastest channel trim)
    % Needs refLbl (same row order as your epochs). Skip if unknown.
    % XTrain2 = keep_grads(XTrain2, refLbl);
    % XVal2   = keep_grads(XVal2,   refLbl);

    % 2) Crop around likely transient (speeds training a lot)
    targetWin = 120;  % samples (adjust to your Fs/window, e.g., ~200 ms)
    XTrain2 = crop_by_peakgrad(XTrain2, targetWin);
    XVal2   = crop_by_peakgrad(XVal2,   targetWin);

    % 3) Low-pass + 2× decimate time (halve T while preserving spike edges)
    Fs   = 1000;    % <-- set to your sampling rate (Hz)
    fcut = 80;      % low-pass cutoff
    XTrain2 = decimate2_lp(XTrain2, Fs, fcut);
    XVal2   = decimate2_lp(XVal2,   Fs, fcut);

    % 4) Spatial PCA to r components (fit on TRAIN only, apply to both)
    r = 64;   % try 4864; often lossless for spikes
    [XTrain2, Wpca] = pca_project_train(XTrain2, r);
    XVal2            = pca_project_apply(XVal2,   Wpca);

    % 5) Fix a single length for both sets (avoid "learn the padding")
    Ttrain = cellfun(@(x) size(x,2), XTrain2);
    Tval   = cellfun(@(x) size(x,2), XVal2);
    Tfix   = max(1, round(prctile([Ttrain(:); Tval(:)], 95)));
    XTrain2 = pad_sequences_to_length(XTrain2, Tfix);
    XVal2   = pad_sequences_to_length(XVal2,   Tfix);

    % (optional) Switch data to single to save memory
    XTrain2 = cellfun(@single, XTrain2, 'UniformOutput', false);
    XVal2   = cellfun(@single, XVal2,   'UniformOutput', false);

    % Update input size for model
    inputSize = size(XTrain2{1},1);  % now r (e.g., 64)
    padArgs   = {"SequenceLength", Tfix, "SequencePaddingDirection","right"};

end

%% ----- Optional: Transformer positional channels (do this BEFORE inputSize)
if strcmpi(opts.arch,'transformer')
    % gentler PE scale (often a bit stabler)
    addPos = @(x) [x; 0.05*make_posenc(size(x,2),32)]; % adds 32 rows
    XTrain2 = cellfun(addPos, XTrain2,'UniformOutput',false);
    XVal2   = cellfun(addPos, XVal2  ,'UniformOutput',false);
end

switch lower(opts.arch)
    case 'resnet1d'
        % ResNet-1D with dilations; good default for MEG IED detection
        mb = 128; lr = 3e-4; maxE = 150;
        valFreq = max(5, ceil(numel(XTrain2)/mb));

        lgraph = resnet1d_layers(inputSize, 2);
        lgraph = replaceLayer(lgraph, "cls", clsLayer);

        trainOpts = trainingOptions('adam', ...
            'InitialLearnRate',lr, ...
            'LearnRateSchedule','piecewise','LearnRateDropPeriod',30,'LearnRateDropFactor',0.5, ...
            'MaxEpochs',maxE, ...
            'MiniBatchSize',mb, ...
            'L2Regularization',1e-3, ...
            'GradientThresholdMethod','l2norm', ...
            'GradientThreshold',1, ...
            padArgs{:}, ...
            'Shuffle','every-epoch', ...
            'ValidationData',{XVal2,YVal}, ...
            'ValidationFrequency',valFreq, ...
            'ValidationPatience',15, ...
            'ExecutionEnvironment', ternary(opts.useGPU,'gpu','cpu'), ...
            'Verbose',opts.verbose, ...
            'Plots','training-progress', ...
            'CheckpointPath',opts.modelDir);

        mdlName = 'resnet1d_spike_classifier.mat';
        net = trainNetwork(XTrain2, YTrain, lgraph, trainOpts);
    case 'resnet_masked'

        [net, dlnet_pre, cm, acc, bal] = run_A3_finetune_masked( ...
    XTrain2, YTrain, XVal2, YVal, ...
    'FixedLenPct', 85, 'MaskFrac', 0.2, ...
    'PreEpochs', 10, 'FT_Epochs', 30, ...
    'MB', 128, 'LR', 3e-4);


        %% 0) (Optional) decimate to cut time in half/quarter
        decim = 2;  % try 2 or 4 if fs is high
        XTrain_fast = cellfun(@(x) single(x(:,1:decim:end)), XTrain2, 'uni', false);
        XVal_fast   = cellfun(@(x) single(x(:,1:decim:end)),   XVal2, 'uni', false);

        %% 1) Tight max sequence length (pads/crops right)
        maxLen = round(prctile(cellfun(@(x) size(x,2), XTrain_fast), 70));  % 6075 is fine

        %% 2) Sizes & class weights
        inputSize = size(XTrain_fast{1},1);
        cats = categories(YTrain);
        YTrain = reordercats(YTrain, cats);
        YVal   = reordercats(YVal,   cats);

        cnt = countcats(YTrain);
        w   = sum(cnt)./max(cnt,1); w = w/mean(w);   % class weights

        %% 3) Tiny sequence network (fast)
        layers = [
            sequenceInputLayer(inputSize, 'Normalization','none','Name','in')
            bilstmLayer(32, 'OutputMode','last','Name','bi1')   % small!
            dropoutLayer(0.2,'Name','drop')
            fullyConnectedLayer(numel(cats),'Name','fc')
            softmaxLayer('Name','sm')
            classificationLayer('Name','cls','Classes',cats,'ClassWeights',w)
            ];

        opts = trainingOptions('sgdm', ...
            'InitialLearnRate',5e-3, 'Momentum',0.9, ...
            'MaxEpochs',20, 'MiniBatchSize',256, ...              % push batch up
            'Shuffle','every-epoch', ...
            'SequenceLength',maxLen, 'SequencePaddingDirection','right', ...
            'GradientThresholdMethod','l2norm','GradientThreshold',1, ...
            'L2Regularization',1e-4, ...
            'ValidationData',{XVal_fast, YVal}, ...
            'ValidationFrequency', max(20, ceil(numel(XTrain_fast)/256)), ...
            'ValidationPatience',5, ...
            'Plots','training-progress', ...
            'ExecutionEnvironment','auto','Verbose',false);

        net = trainNetwork(XTrain_fast, YTrain, layers, opts);

    case 'bilstm_with_regularization'

        % ----- Model (smaller + more dropout) -----
        layers = [
            sequenceInputLayer(inputSize,'Normalization','none','Name','in')
            bilstmLayer(64,'OutputMode','last','Name','bilstm1')   % was 128/64 & 2 layers
            dropoutLayer(0.35,'Name','drop1')                      % ? dropout
            fullyConnectedLayer(32,'Name','fc1')
            reluLayer('Name','relu1')
            dropoutLayer(0.35,'Name','drop2')
            fullyConnectedLayer(2,'Name','fc2')
            softmaxLayer('Name','sm')
            classificationLayer('Name','cls')];

        % ----- Training options (lower LR, weight decay, LR schedule, tighter ES) -----
        mb      = 64;
        valFreq = max(10, ceil(numel(XTrain2)/mb));
        lr0     = 3e-4;

        trainOpts = trainingOptions('adam', ...
            'InitialLearnRate', lr0, ...
            'LearnRateSchedule','piecewise','LearnRateDropPeriod',25,'LearnRateDropFactor',0.5, ...
            'L2Regularization',5e-4, ...              % weight decay
            'MaxEpochs',120, ...
            'MiniBatchSize',mb, ...
            'GradientThresholdMethod','l2norm','GradientThreshold',1, ...
            'Shuffle','every-epoch', ...
            padArgs{:}, ...                           % use fixed SequenceLength (not "longest")
            'ValidationData',{XVal2,YVal}, ...
            'ValidationFrequency',valFreq, ...
            'ValidationPatience',10, ...              % earlier stop
            'ExecutionEnvironment', ternary(opts.useGPU,'gpu','cpu'), ...
            'Verbose',opts.verbose, ...
            'Plots','training-progress', ...
            'CheckpointPath',opts.modelDir);

        % If you computed class weights -> replace the cls layer before training:
        % idx = arrayfun(@(L) isa(L,'nnet.cnn.layer.ClassificationLayer'), layers);
        % layers(idx) = clsLayer;
        % swap in weighted classifier
        idx = arrayfun(@(L) isa(L,'nnet.cnn.layer.ClassificationLayer'), layers);
        layers(idx) = clsLayer;

        mdlName = 'bilstm_reg_spike_classifier.mat';


        net = trainNetwork(XTrain2, YTrain, layers, trainOpts);
    case 'longview'

        % Build fast long-view features (vectorized, no parfor)
        XTrain7 = cell(size(XTrain2));
        for i = 1:numel(XTrain2)
            XTrain7{i} = fast_longview_features(single(XTrain2{i}), fs);  % -> [C x T x D]
        end
        XVal7 = cell(size(XVal2));
        for i = 1:numel(XVal2)
            XVal7{i} = fast_longview_features(single(XVal2{i}), fs);
        end

        % Pack to 4-D arrays with a fixed width for both train/val
        fixedLen = round(prctile(cellfun(@(x) size(x,2), XTrain7), 90));
        [trainImgs, trainLabs] = cell2imgs(XTrain7, YTrain, fixedLen);  % [C x W x D x N]
        [valImgs,   valLabs]   = cell2imgs(XVal7,   YVal,   fixedLen);

        % Sizes and depth D come directly from the packed arrays
        C = size(trainImgs,1);
        W = size(trainImgs,2);
        D = size(trainImgs,3);   % <-- use this; no repmat

        % Class weights (pre-oversample counts)
        cats = categories(trainLabs);
        cnt  = countcats(trainLabs);
        w    = sum(cnt)./max(cnt,1); w = w/mean(w);
        clsLayer = classificationLayer('Name','cls','Classes',cats,'ClassWeights',w);

        % Build LV-CadeNet-Lite graph (uses D from above)
        lg = lv_cadenet_lite_layers(C, W, D, 2, clsLayer);

        % Train
        mb = 128; valFreq = max(10, ceil(size(trainImgs,4)/mb));
        opts = trainingOptions('adam', ...
            'InitialLearnRate',3e-4, 'MaxEpochs',60, 'MiniBatchSize',mb, ...
            'Shuffle','every-epoch', ...
            'ValidationData',{valImgs, valLabs}, ...
            'ValidationFrequency', valFreq, 'ValidationPatience',10, ...
            'L2Regularization',1e-4, ...
            'ExecutionEnvironment','auto', ...
            'Plots','training-progress', 'Verbose',true);

        net = trainNetwork(trainImgs, trainLabs, lg, opts);

    case 'lstm'
        mb = 64; lr = 3e-4; maxE = 120;              % lower LR is gentler
        valFreq = max(5, ceil(numel(XTrain2)/mb));

        layers = [
            sequenceInputLayer(inputSize,'Normalization','none')
            bilstmLayer(256,'OutputMode','sequence')
            dropoutLayer(0.2)
            bilstmLayer(64,'OutputMode','last')
            dropoutLayer(0.2)
            fullyConnectedLayer(2)
            softmaxLayer
            classificationLayer('Name','cls')];      % <-- name it here

        % Swap in weighted classifier
        idx = arrayfun(@(L) isa(L,'nnet.cnn.layer.ClassificationLayer'), layers);
        layers(idx) = clsLayer;

        trainOpts = trainingOptions('adam', ...
            'InitialLearnRate', lr, ...
            'LearnRateSchedule','piecewise','LearnRateDropPeriod',25,'LearnRateDropFactor',0.5, ... % optional
            'L2Regularization',5e-4, ...                          % optional weight decay
            'MaxEpochs', maxE, ...
            'MiniBatchSize', mb, ...
            'GradientThresholdMethod','l2norm','GradientThreshold',1, ...
            'Shuffle','every-epoch', ...
            padArgs{:}, ...
            'ValidationData',{XVal2,YVal}, ...
            'ValidationFrequency', valFreq, ...
            'ValidationPatience', 10, ...
            'ExecutionEnvironment', ternary(opts.useGPU,'gpu','cpu'), ...
            'Verbose', opts.verbose, ...
            'Plots','training-progress', ...
            'CheckpointPath', opts.modelDir);

        net = trainNetwork(XTrain2, YTrain, layers, trainOpts);
        mdlName = 'lstm_spike_classifier.mat';

    case 'gru_td' % GRU with time-distributed dense via fold/unfold
        mb = 64; lr = 1e-3; maxE = 100;
        valFreq = max(5, ceil(numel(XTrain2)/mb));
        lg = layerGraph();
        lg = addLayers(lg, sequenceInputLayer(inputSize,"Normalization","none","Name","input"));
        lg = addLayers(lg, sequenceFoldingLayer("Name","fold"));
        lg = addLayers(lg, fullyConnectedLayer(32,"Name","fc_td"));
        lg = addLayers(lg, reluLayer("Name","relu_td"));
        lg = addLayers(lg, sequenceUnfoldingLayer("Name","unfold"));
        lg = addLayers(lg, flattenLayer("Name","flat"));
        lg = addLayers(lg, [
            gruLayer(64,"OutputMode","last","Name","gru")
            dropoutLayer(0.2,"Name","drop")
            fullyConnectedLayer(2,"Name","fc")
            softmaxLayer("Name","sm")
            classificationLayer("Name","cls")]);
        lg = connectLayers(lg,"input","fold/in");
        lg = connectLayers(lg,"fold/out","fc_td");
        lg = connectLayers(lg,"fc_td","relu_td");
        lg = connectLayers(lg,"relu_td","unfold/in");
        lg = connectLayers(lg,"unfold/out","flat");
        lg = connectLayers(lg,"flat","gru");
        lg = connectLayers(lg,"fold/miniBatchSize","unfold/miniBatchSize");
        trainOpts = trainingOptions("adam", ...
            "InitialLearnRate",lr, ...
            "MaxEpochs",maxE, ...
            "MiniBatchSize",mb, ...
            'GradientThresholdMethod','l2norm', ...
            'GradientThreshold',1, ...
            padArgs{:}, ...
            "ValidationData",{XVal2,YVal}, ...
            "ValidationFrequency",valFreq, ...
            "ValidationPatience",15, ...
            "Shuffle","every-epoch", ...
            "ExecutionEnvironment", ternary(opts.useGPU,'gpu','cpu'), ...
            "Verbose",opts.verbose, ...
            "Plots","training-progress", ...
            "CheckpointPath",opts.modelDir);
        lg = replaceLayer(lg, "cls", clsLayer);
        mdlName = 'grutd_spike_classifier.mat';
        net = trainNetwork(XTrain2,YTrain,lg,trainOpts);

    case 'transformer'
        mb = 128; lr = 3e-4; maxE = 150;
        valFreq = max(5, ceil(numel(XTrain2)/mb));
        dModel  = opts.dModel; numHeads = opts.numHeads; keyDim = dModel/numHeads; % integer by assert in parse_inputs
        lg = layerGraph();
        lg = addLayers(lg, sequenceInputLayer(inputSize,"Normalization","none","Name","in"));
        lg = addLayers(lg, fullyConnectedLayer(dModel,"Name","proj"));
        % Encoder (pre-norm)
        lg = addLayers(lg, layerNormalizationLayer("Name","ln1"));
        lg = addLayers(lg, selfAttentionLayer(numHeads,keyDim, "DropoutProbability",0.1,"Name","mha1"));
        lg = addLayers(lg, dropoutLayer(0.1,"Name","drop_attn"));
        lg = addLayers(lg, additionLayer(2,"Name","add_attn"));
        lg = addLayers(lg, layerNormalizationLayer("Name","ln2"));
        lg = addLayers(lg, [
            fullyConnectedLayer(4*dModel,"Name","ff1")
            reluLayer("Name","relu1")
            dropoutLayer(0.1,"Name","drop_ff")
            fullyConnectedLayer(dModel,"Name","ff2")]);
        lg = addLayers(lg, additionLayer(2,"Name","add_ff"));
        % Head
        lg = addLayers(lg, globalAveragePooling1dLayer("Name","gap"));
        lg = addLayers(lg, dropoutLayer(0.2,"Name","head_drop"));
        lg = addLayers(lg, fullyConnectedLayer(2,"Name","fc"));
        lg = addLayers(lg, softmaxLayer("Name","sm"));
        lg = addLayers(lg, classificationLayer("Name","cls"));
        % wires
        lg = connectLayers(lg,"in","proj");
        lg = connectLayers(lg,"proj","ln1");
        lg = connectLayers(lg,"ln1","mha1");
        lg = connectLayers(lg,"mha1","drop_attn");
        lg = connectLayers(lg,"drop_attn","add_attn/in1");
        lg = connectLayers(lg,"ln1","add_attn/in2");
        lg = connectLayers(lg,"add_attn","ln2");
        lg = connectLayers(lg,"ln2","ff1");
        lg = connectLayers(lg,"ff2","add_ff/in1");
        lg = connectLayers(lg,"add_attn","add_ff/in2");
        lg = connectLayers(lg,"add_ff","gap");
        lg = connectLayers(lg,"gap","head_drop");
        lg = connectLayers(lg,"head_drop","fc");
        lg = connectLayers(lg,"fc","sm");
        lg = connectLayers(lg,"sm","cls");
        trainOpts = trainingOptions("adam", ...
            "InitialLearnRate",lr, ...
            "MaxEpochs",maxE, ...
            "MiniBatchSize",mb, ...
            'GradientThresholdMethod','l2norm', ...
            'GradientThreshold',1, ...
            padArgs{:}, ...
            "Shuffle","every-epoch", ...
            "ValidationData",{XVal2,YVal}, ...
            "ValidationFrequency",valFreq, ...
            "ValidationPatience",15, ...
            "ExecutionEnvironment", ternary(opts.useGPU,'gpu','cpu'), ...
            "Verbose",opts.verbose, ...
            "Plots","training-progress", ...
            "CheckpointPath",opts.modelDir);
        lg = replaceLayer(lg, "cls", clsLayer);
        net = trainNetwork(XTrain2,YTrain,lg,trainOpts);
        mdlName = 'transformer_spike_classifier.mat';

    case 'ems'

        fixedLen = round(prctile(cellfun(@(x) size(x,2), XTrain2), 90));   % tighter than 95
        [trainImgs, trainLabs] = cell2imgs(XTrain2, YTrain, fixedLen);
        [valImgs,   valLabs]   = cell2imgs(XVal2,   YVal,   fixedLen);

        megH = size(trainImgs,1); megW = size(trainImgs,2);
        lg = emsnet_lite_layers(megH, megW, 2, clsLayer);

        mb = 128;
        trainOpts = trainingOptions('adam', ...
            'InitialLearnRate',3e-4, 'MaxEpochs',60, 'MiniBatchSize',mb, ...
            'Shuffle','every-epoch', ...
            'ValidationData',{valImgs, valLabs}, ...
            'ValidationFrequency', max(10, ceil(size(trainImgs,4)/mb)), ...
            'ValidationPatience',10, ...
            'ExecutionEnvironment','auto', ...  % GPU if available
            'Plots','training-progress', ...
            'Verbose',true);
        net = trainNetwork(trainImgs, trainLabs, lg, trainOpts);

    otherwise
        error('Unknown arch: %s. Use ''lstm'', ''gru_td'', or ''transformer''.',opts.arch);
end

%% ------------------------------ SAVE & EVAL -----------------------------
% Predictions + metrics
[YP, score] = classify(net,XVal2);
cm = confusionmat(YVal, YP,'Order',categories(YVal));
acc = sum(diag(cm))/sum(cm,'all');
prec = diag(cm)./max(1,sum(cm,1)');
rec  = diag(cm)./max(1,sum(cm,2));
f1   = 2*(prec.*rec)./max(prec+rec,eps);
balAcc = mean(rec,'omitnan');

% ROC (optional if Statistics Toolbox present)
auc = NaN; rocX = []; rocY = [];
posClass = 'Spike';
if exist('perfcurve','file')==2
    [~,posIdx] = ismember(posClass, categories(YVal));
    scoresPos = score(:,posIdx);
    [rocX,rocY,~,auc] = perfcurve(YVal, scoresPos, posClass);
end

% Save reproducibility meta (git hash optional)
[st,hash] = system(sprintf('git -C "%s" rev-parse HEAD', opts.codeRoot));
if st==0, git_hash = strtrim(hash); else, git_hash='NA'; end
% meta = struct('matlabVersion',version,'gitHash',git_hash,'timestamp',string(datetime("now"),'yyyy-mm-dd HH:MM:SS'));

meta = struct('matlabVersion',version, ...
              'gitHash',git_hash, ...
              'timestamp', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));

save(fullfile(opts.modelDir,mdlName), 'net','normStats','mu','sd','opts','trainOpts','cm','acc','balAcc','f1','auc','meta');
save(fullfile(opts.modelDir,'run_meta.mat'),'meta');

% Figures
if opts.savePlots
    fig1=figure; plotconfusion(YVal,YP); title(sprintf('Val Acc=%.2f, BA=%.2f, AUC=%.2f',acc,balAcc,auc));
    saveas(fig1, fullfile(opts.modelDir,'confusion_val.png'));
    % close(fig1);
    if ~isempty(rocX)
        fig2=figure; plot(rocX,rocY); xlabel('FPR'); ylabel('TPR'); title(sprintf('ROC (AUC=%.3f)',auc)); grid on;
        saveas(fig2, fullfile(opts.modelDir,'roc_val.png'));
        % close(fig2);
    end
    if ~isempty(XTrain2)
        fig3=figure('Name','Example epoch'); clf
        plot(XTrain2{1}'); title('First normalized epoch'); xlabel('Time');
        saveas(fig3, fullfile(opts.modelDir,'example_epoch.png'));
        % close(fig3);
    end
end

end % main

% ========================= HELPER FUNCTIONS =============================
function [cells,groups,ref] = build_cells(files,doInit,groupFcn)
% Converts anot_data_all entries into 2-D epoch matrices; collects group IDs.
if nargin<2, doInit = 0; end
if nargin<3, groupFcn = @extract_group_id; end

cells = {}; groups = {}; ref = [];
totalFiles = numel(files); totEpochs = 0;

ft_progress('init','text', 'Reading %d files...', totalFiles);
c = onCleanup(@() ft_progress('close')); % ensure close even on error
for k = 1:totalFiles
    S = load(files(k).name);
    if ~isfield(S,'anot_data_all'), continue; end
    groupID = groupFcn(files(k).name);
    for t = 1:numel(S.anot_data_all)
        D = S.anot_data_all{t};
        if doInit && isempty(ref), ref = D.label; end
        D = do_ensure_consistent_labels(D, ref);
        allLbl  = D.label;
        megMask = startsWith(allLbl,'MEG');   % drop EEG, misc
        MEG_labels = allLbl(megMask);
        E  = D.trial{1}(megMask ,:);          % (channels x time)
        E  = maggrad_scale(E,MEG_labels);
        cells{end+1,1}  = do_normalize_data(E,'demean'); %#ok<AGROW>
        groups{end+1,1} = groupID;                       %#ok<AGROW>
        totEpochs = totEpochs + 1;
    end
    ft_progress(k/totalFiles, 'File %d/%d | epochs read: %d', k, totalFiles, totEpochs);
end
fprintf('Finished build_cells: %d files, %d total epochs.\n', totalFiles, totEpochs);
end

function [cells,groups] = build_cells_with_ref(files,ref,groupFcn)
% Build cells aligned to a provided reference label list
if nargin<3, groupFcn = @extract_group_id; end
cells = {}; groups = {};
totalFiles = numel(files); totEpochs = 0;
ft_progress('init','text','Reading %d files (aligned to ref)...', totalFiles);
c = onCleanup(@() ft_progress('close'));
for k=1:totalFiles
    S = load(files(k).name);
    if ~isfield(S,'anot_data_all'), continue; end
    groupID = groupFcn(files(k).name);
    for t=1:numel(S.anot_data_all)
        D = do_ensure_consistent_labels(S.anot_data_all{t}, ref);
        allLbl  = D.label;
        megMask = startsWith(allLbl,'MEG');
        MEG_labels = allLbl(megMask);
        E = D.trial{1}(megMask,:);
        E = maggrad_scale(E,MEG_labels);
        cells{end+1,1} = do_normalize_data(E,'demean'); %#ok<AGROW>
        groups{end+1,1} = groupID; %#ok<AGROW>
        totEpochs = totEpochs + 1;
    end
    ft_progress(k/totalFiles, 'File %d/%d | epochs read: %d', k, totalFiles, totEpochs);
end
fprintf('Finished build_cells_with_ref: %d files, %d total epochs.', totalFiles, totEpochs);
end

function files = list_mat_files(root)
% Robust file discovery without rdir dependency
if exist('rdir','file')==2
    files = rdir(fullfile(root,'*.mat'));
else
    dd = dir(fullfile(root,'**','*.mat')); % requires R2016b+
    for i=1:numel(dd), dd(i).name = fullfile(dd(i).folder, dd(i).name); end
    files = dd;
end
end

function gid = extract_group_id(fullpath)
% Heuristics to produce a subject/session group ID from the path/filename
[folder, base, ~] = fileparts(fullpath);
[~, parent] = fileparts(folder);
m = regexp(base,'^([^_]+)_','tokens','once');
if ~isempty(m), gid = m{1}; else, gid = parent; end
end

% function [trainMask,valMask] = group_holdout_split(groups, valFrac, seed)
% % Splits by unique group IDs to prevent leakage between train/val
% if nargin<3, seed = 0; end
% rng(seed,'twister');
% u = unique(groups,'stable');
% G = numel(u);
% if G == 0
%     error('group_holdout_split:NoGroups','No groups provided.');
% elseif G == 1
%     % Put the single group entirely into TRAIN; keep VAL empty
%     valMask = false(size(groups));
%     trainMask = true(size(groups));
%     return;
% end
% nVal = round(valFrac * G);
% % Ensure at least 1 val group and at least 1 train group
% nVal = max(1, min(nVal, G-1));
% idx = randperm(G);
% valGroups = u(idx(1:nVal));
% valMask = ismember(groups, valGroups);
% trainMask = ~valMask;
% % end
% rng(seed,'twister');
% u = unique(groups,'stable');
% nVal = max(1, round(valFrac * numel(u)));
% idx = randperm(numel(u));
% valGroups = u(idx(1:nVal));
% valMask = ismember(groups, valGroups);
% trainMask = ~valMask;
% end

function [trainMask,valMask] = group_holdout_split(groups, valFrac, seed)
if nargin<3, seed = 0; end
rng(seed,'twister');
u = unique(groups,'stable'); G = numel(u);
if G==0, error('group_holdout_split:NoGroups','No groups provided.'); end
if G==1, valMask = false(size(groups)); trainMask = true(size(groups)); return; end
nVal = max(1, min(round(valFrac*G), G-1));
idx = randperm(G);
valGroups = u(idx(1:nVal));
valMask   = ismember(groups, valGroups);
trainMask = ~valMask;
end


function [mu,sd] = fit_channel_stats(C)
% Compute per-channel mean & std over variable-length sequences
ch = size(C{1},1); sumC = zeros(ch,1); sumSqC = zeros(ch,1); nTot=0;
for i=1:numel(C)
    xi=double(C{i}); sumC=sumC+sum(xi,2); sumSqC=sumSqC+sum(xi.^2,2); nTot=nTot+size(xi,2);
end
mu = sumC./nTot; varC = (sumSqC./nTot) - mu.^2; sd = sqrt(max(varC,1e-12));
end

function [med,madv] = fit_channel_stats_robust(C)
if isempty(C) || isempty(C{1})
    error('fit_channel_stats_robust:EmptyTraining','XTrain is empty; cannot fit robust stats.');
end
ch = size(C{1},1); med = zeros(ch,1); madv = ones(ch,1);
for r = 1:ch
    rows = cellfun(@(z) double(z(r,:)), C, 'UniformOutput', false);
    rows = rows(~cellfun(@isempty, rows));
    if isempty(rows), med(r)=0; madv(r)=1; continue; end
    x = [rows{:}];
    if isempty(x) || all(isnan(x))
        med(r)=0; madv(r)=1;
    else
        m = median(x,'omitnan');
        med(r)=m;
        madv(r)=max(1e-12, 1.4826*median(abs(x-m),'omitnan'));
    end
end
end

% function [med,madv] = fit_channel_stats_robust(C)
% % Robust per-channel stats (median, MADsigma) for variable-length epochs
% % Avoid cell2mat on unequal lengths; concatenate row-vectors horizontally.
% if isempty(C) || isempty(C{1})
%     error('fit_channel_stats_robust:EmptyTraining','XTrain is empty; cannot fit robust stats.');
% end
% ch = size(C{1},1); med = zeros(ch,1); madv = zeros(ch,1);
% for r = 1:ch
%     rows = cellfun(@(z) double(z(r,:)), C, 'UniformOutput', false);
%     rows = cellfun(@(v) v(:).', rows, 'UniformOutput', false); % ensure row vectors
%     rows = rows(~cellfun(@isempty, rows));
%     if isempty(rows)
%         med(r)  = 0;      % fallback
%         madv(r) = 1;      % fallback
%         continue;
%     end
%     x = [rows{:}];       % 1 x (sum T)
%     if isempty(x) || all(isnan(x))
%         med(r)  = 0; madv(r) = 1;
%     else
%         m = median(x,'omitnan');
%         med(r)  = m;
%         madv(r) = 1.4826 * median(abs(x - m),'omitnan');
%     end
% end
% madv = max(madv,1e-12);  % protect against divide-by-zero
% % end
% x = [rows{:}];       % 1 x (sum T) vector
% if isempty(x) || all(isnan(x))
%     med(r)  = 0; madv(r) = 1;  % safe defaults
% else
%     m = median(x,'omitnan');
%     med(r)  = m;
%     madv(r) = 1.4826 * median(abs(x - m),'omitnan');  % MAD -> ~sigma
% end
% % end
% % madv = max(madv,1e-12);  % protect against divide-by-zero
% % end
% madv = max(madv,1e-12);
% end

function [Xb,Yb] = balance_oversample(X,Y)
% Oversample minority class (simple, version-proof)
cats = categories(Y); counts = countcats(Y);
[~,maj] = max(counts); [~,minc] = min(counts);
majLab = cats{maj}; minLab = cats{minc};
minIdx = find(Y==minLab); majIdx = find(Y==majLab);
rep = numel(majIdx) - numel(minIdx);
if rep <= 0, Xb = X; Yb = Y; return; end
addIdx = minIdx(randi(numel(minIdx), rep, 1));     % sample w/ replacement
Xb = [X; X(addIdx)];  Yb = [Y; Y(addIdx)];
perm = randperm(numel(Yb)); Xb = Xb(perm); Yb = Yb(perm);
end

function aug = default_meg_aug()
% Default augmentation hyperparameters (replicates the paper's four augs)
aug = struct();
aug.p_crop   = 0.8;     % probability to random-crop
aug.crop_len = [];      % if empty, random fraction in [0.5, 1.0]
aug.crop_frac_min = 0.5;
aug.crop_frac_max = 1.0;
aug.p_shift  = 0.5;     % probability to sensor-order shift
aug.max_shift= 16;      % max circular shift (channels)
aug.p_flip   = 0.5;     % probability to phase flip (× -1)
aug.p_mix    = 0.5;     % probability to do frequency mix with a NoSpike partner
aug.mix_alpha_min = 0.15;  % how much partner magnitude to inject
aug.mix_alpha_max = 0.35;
end

function aug = fill_default_meg_aug(user)
% Merge user-provided overrides with defaults
aug = default_meg_aug();
uf = fieldnames(user);
for i=1:numel(uf), aug.(uf{i}) = user.(uf{i}); end
end

function A = make_meg_augmentor(aug, negPool, chanInfo)
% Returns a function handle A(x) that applies the four augmentations
if nargin<3 || isempty(chanInfo), chanInfo = struct(); end
if nargin<2 || isempty(negPool),  negPool  = {}; end

A = @augment_one;

    function y = augment_one(x)
        y = x;
        % 1) Cropping window (random)
        if rand < aug.p_crop
            y = aug_crop_window(y, aug);
        end
        % 2) Sensor-order shift (label-aware if masks provided)
        if rand < aug.p_shift
            y = aug_channel_shift(y, chanInfo, aug.max_shift);
        end
        % 3) Phase flip
        if rand < aug.p_flip
            y = -y;
        end
        % 4) Frequency mixing with a negative partner (label-preserving)
        if ~isempty(negPool) && rand < aug.p_mix
            partner = negPool{randi(numel(negPool))};
            y = aug_mix_frequency(y, partner, aug);
        end
    end
end

function y = aug_crop_window(x, aug)
% Random crop along time (keeps variable length; padding happens in training)
T = size(x,2);
if isempty(aug.crop_len)
    frac = aug.crop_frac_min + (aug.crop_frac_max-aug.crop_frac_min)*rand;
    w = max(1, min(T, round(frac*T)));
else
    w = max(1, min(T, round(aug.crop_len)));
end
s = randi([1, max(1, T-w+1)]);
y = x(:, s:s+w-1);
end

function y = aug_channel_shift(x, chanInfo, maxShift)
% Circularly shift channel order; if mag/grad masks available, shift within each set
if nargin<3 || isempty(maxShift), maxShift = 16; end
n = size(x,1); y = x;
if isfield(chanInfo,'idxMag') && ~isempty(chanInfo.idxMag)
    idx = chanInfo.idxMag(:)'; k = randi([-maxShift,maxShift]);
    order = circshift(1:numel(idx), [0 k]);
    y(idx,:) = x(idx(order),:);
end
if isfield(chanInfo,'idxGrad') && ~isempty(chanInfo.idxGrad)
    idx = chanInfo.idxGrad(:)'; k = randi([-maxShift,maxShift]);
    order = circshift(1:numel(idx), [0 k]);
    y(idx,:) = x(idx(order),:);
end
if (~isfield(chanInfo,'idxMag') || isempty(chanInfo.idxMag)) && ...
        (~isfield(chanInfo,'idxGrad') || isempty(chanInfo.idxGrad))
    k = randi([-maxShift,maxShift]);
    order = circshift(1:n, [0 k]);
    y = x(order,:);
end
end

function y = aug_mix_frequency(x, partner, aug)
% Mix magnitude spectra of x with a partner (keep x phase); channelwise
T = min(size(x,2), size(partner,2));
if T <= 1
    y = x; return;
end
x = x(:,1:T); p = partner(:,1:T);
X = fft(x, [], 2); P = fft(p, [], 2);
alpha = aug.mix_alpha_min + (aug.mix_alpha_max-aug.mix_alpha_min)*rand;
mag = (1-alpha)*abs(X) + alpha*abs(P);
phs = angle(X);
Z = mag .* exp(1j*phs);
y = real(ifft(Z, [], 2));
% end
% Xaug = x;
end

function v = ensure_column(v)
% Force a column vector shape
v = v(:);
end

function X = maggrad_scale(X,chLabels)
% Scale magnetometers and gradiometers to comparable ranges
isMag  = endsWith(chLabels,'1');
isGrad = endsWith(chLabels,'2') | endsWith(chLabels,'3');
mMag  = median(abs(X(isMag ,:)), 'all','omitnan') + eps;
mGrad = median(abs(X(isGrad,:)), 'all','omitnan') + eps;
X(isMag ,:) = X(isMag ,:) ./ mMag;
X(isGrad,:) = X(isGrad,:) ./ mGrad;
end

function PE = make_posenc(T,d)
% Sin/cos positional channels, returned as [d x T], d must be even
if mod(d,2)~=0, error('d must be even'); end
t = (0:T-1)';                        % T x 1
i = (0:(d/2-1));                    % 1 x (d/2)
angles = t ./ (10000.^(2*i/d));     % T x (d/2)
PE = [sin(angles) cos(angles)]';    % d x T
end

function out = parse_inputs(def, varargin)
% Robust parser: accepts 'norm','augment','aug','seed', etc.
p = inputParser;
p.FunctionName   = mfilename;          % show the real function name in errors
p.CaseSensitive  = false;
p.PartialMatching= true;
p.KeepUnmatched  = false;

% Existing params
addParameter(p,'ftRoot',   def.ftRoot);
addParameter(p,'codeRoot', def.codeRoot);
addParameter(p,'spikeDir', def.spikeDir);
addParameter(p,'noSpikeDir',def.noSpikeDir);
addParameter(p,'skipList', def.skipList);
addParameter(p,'useGPU',   def.useGPU);
addParameter(p,'modelDir', def.modelDir);
addParameter(p,'savePlots',def.savePlots);
addParameter(p,'verbose',  def.verbose);
addParameter(p,'valFrac',  def.valFrac);

% Model selection
% validArch = @(s) any(strcmpi(s,{'lstm','gru_td','transformer'}));
validArch = @(s) any(strcmpi(s,{'lstm','gru_td','transformer','resnet1d'}));  % < added 'resnet1d'
addParameter(p,'arch',    def.arch,    validArch);
addParameter(p,'dModel',  def.dModel,  @(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
addParameter(p,'numHeads',def.numHeads,@(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));

% Normalization + augmentation
addParameter(p,'norm',    def.norm,    @(s) any(strcmpi(s,{'zscore','robust'})));
addParameter(p,'augment', def.augment, @(x) (islogical(x)||isnumeric(x)) && isscalar(x));
addParameter(p,'aug',     struct(),    @(s) isempty(s) || isstruct(s));   % <-- NEW: accept aug struct

% Repro seed
addParameter(p,'seed',    def.seed,    @(x) isnumeric(x) && isscalar(x) && isfinite(x));

parse(p, varargin{:});
out = p.Results;

% Post-parse sanitation
out.augment = logical(out.augment);  % accept 0/1 or true/false

% Transformer integrity check
if strcmpi(out.arch,'transformer')
    assert(mod(out.dModel,out.numHeads)==0, ...
        'dModel (%d) must be divisible by numHeads (%d).', out.dModel, out.numHeads);
end
% Do NOT call rng here; the caller sets RNG once via set_rng(out.seed)
end


% function setup_environment(o)
% addpath(o.ftRoot); ft_defaults;
% addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/functions_new')));
% addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/helper')));
% end

function setup_environment(o)
addpath(o.ftRoot); ft_defaults;
addpath(genpath(fullfile(o.codeRoot,'FT_functions','functions_new')));
addpath(genpath(fullfile(o.codeRoot,'FT_functions','helper')));
end

function r = ternary(c,a,b); if c, r=a; else, r=b; end; end

function lgraph = resnet1d_layers(inputSize, numClasses)
% Build a compact ResNet-1D with dilated convolutions for sequence data.
% Input is a sequence: [channels x time], provided as cell array elements.

% Stem
lgraph = layerGraph();
lgraph = addLayers(lgraph, sequenceInputLayer(inputSize, 'Normalization','none', 'Name','in'));
stem = [
    convolution1dLayer(7, 32, 'Padding','same', 'Stride',1, 'Name','stem_conv', 'WeightsInitializer','he')
    batchNormalizationLayer('Name','stem_bn')
    reluLayer('Name','stem_relu')
    ];
lgraph = addLayers(lgraph, stem);
lgraph = connectLayers(lgraph, 'in', 'stem_conv');

prevName = 'stem_relu'; prevCh = 32;

% Residual block schedule (filters, dilation). Nine blocks total.
cfg = [ ...
    64  1
    64  2
    64  4
    128 1
    128 2
    128 4
    256 1
    256 2
    256 4];

for i = 1:size(cfg,1)
    filters = cfg(i,1);
    dil     = cfg(i,2);
    blkName = sprintf('b%d',i);
    [lgraph, prevName, prevCh] = addResBlock1d(lgraph, prevName, prevCh, filters, 3, dil, blkName, 0.1);
end

% Head
head = [
    globalAveragePooling1dLayer("Name","gap")
    dropoutLayer(0.2,"Name","head_drop")
    fullyConnectedLayer(numClasses,"Name","fc")
    softmaxLayer("Name","sm")
    classificationLayer("Name","cls")];
lgraph = addLayers(lgraph, head);
lgraph = connectLayers(lgraph, prevName, 'gap');
end


function [lgraph, outName, outCh] = addResBlock1d(lgraph, inName, inCh, filters, k, dil, blk, dropP)
% Residual block:
% main: Conv1D(k,filters,dil) -> BN -> ReLU -> Dropout -> Conv1D(k,filters) -> BN
% skip: identity or 1x1 Conv1D -> BN
% out : Add(main, skip) -> ReLU

% --- main path (auto-connected inside this array) ---
main = [
    convolution1dLayer(k, filters, 'Padding','same', 'DilationFactor',dil, ...
    'Name',[blk '_conv1'], 'WeightsInitializer','he')
    batchNormalizationLayer('Name',[blk '_bn1'])
    reluLayer('Name',[blk '_relu1'])
    dropoutLayer(dropP,'Name',[blk '_drop'])
    convolution1dLayer(k, filters, 'Padding','same', 'DilationFactor',1, ...
    'Name',[blk '_conv2'], 'WeightsInitializer','he')
    batchNormalizationLayer('Name',[blk '_bn2'])];
lgraph = addLayers(lgraph, main);

% --- skip path (projection only if channels change) ---
if inCh ~= filters
    skip = [
        convolution1dLayer(1, filters, 'Padding','same', ...
        'Name',[blk '_proj'], 'WeightsInitializer','he')
        batchNormalizationLayer('Name',[blk '_bnskip'])];
    lgraph = addLayers(lgraph, skip);
    skipOut = [blk '_bnskip'];
    % only connect input to the head of the skip path
    lgraph = connectLayers(lgraph, inName, [blk '_proj']);
else
    skipOut = inName; % identity skip
end

% --- addition + output ReLU ---
addL  = additionLayer(2,'Name',[blk '_add']);
relu2 = reluLayer('Name',[blk '_relu2']);
lgraph = addLayers(lgraph, addL);
lgraph = addLayers(lgraph, relu2);

% Wire entry to main head (internals are auto-wired by addLayers)
lgraph = connectLayers(lgraph, inName,        [blk '_conv1']);
% Wire main tail to add/in1
lgraph = connectLayers(lgraph, [blk '_bn2'],  [blk '_add/in1']);
% Wire skip output (proj or identity) to add/in2
lgraph = connectLayers(lgraph, skipOut,       [blk '_add/in2']);
% Wire add to output ReLU
lgraph = connectLayers(lgraph, [blk '_add'],  [blk '_relu2']);

outName = [blk '_relu2'];
outCh   = filters;
end

function [Xtr_s, Xva_s] = smooth_train_val(Xtr, Xva, Fs, cfg)
Xtr_s = cellfun(@(x) smooth_epoch(x, Fs, cfg), Xtr, 'UniformOutput', false);
Xva_s = cellfun(@(x) smooth_epoch(x, Fs, cfg), Xva, 'UniformOutput', false);
end

function Xs = smooth_epoch(X, Fs, cfg)
Xs = X; if isempty(X), return; end
hasFs = ~isempty(Fs) && isnumeric(Fs) && isfinite(Fs);

switch lower(cfg.type)
    case 'movmean'
        w = getWin(cfg, hasFs, Fs, 'win_ms','win_samp');
        Xs = movmean(X, w, 2, 'Endpoints','shrink');

    case 'median'
        w = getWin(cfg, hasFs, Fs, 'win_ms','win_samp');
        if mod(w,2)==0, w = w+1; end
        Xs = medfilt1(X.', w, [], 1, 'truncate').';

    case 'sgolay'
        w = getWin(cfg, hasFs, Fs, 'frame_ms','frame_samp');
        if mod(w,2)==0, w = w+1; end
        ord = 3;
        if isfield(cfg,'order') && ~isempty(cfg.order), ord = cfg.order; end
        ord = min(ord, w-1);
        Xs  = sgolayfilt(X.', ord, w, [], 1).';

    case 'gauss'
        sig = getSigma(cfg, hasFs, Fs);                 % samples
        W   = getWin(cfg, hasFs, Fs, 'win_ms','win_samp');
        if isempty(W), W = max(5, 2*ceil(4*sig)+1); end % ~±4s
        if mod(W,2)==0, W = W+1; end
        n = -(W-1)/2:(W-1)/2;
        k = exp(-(n.^2)/(2*sig^2)); k = k/sum(k);
        Xs = conv2(X, k, 'same');

    case 'lowpass'
        assert(hasFs, 'lowpass with fc (Hz) requires Fs. Use type="lowpass_norm" instead.');
        fc  = cfg.fc;
        ord = 4; if isfield(cfg,'order') && ~isempty(cfg.order), ord = cfg.order; end
        [b,a] = butter(ord, fc/(Fs/2), 'low');
        Xs = filtfilt(b, a, double(X.')).';

    case 'lowpass_norm'
        % normalized cutoff Wn in (0,1)
        Wn  = cfg.Wn;
        ord = 4; if isfield(cfg,'order') && ~isempty(cfg.order), ord = cfg.order; end
        [b,a] = butter(ord, Wn, 'low');
        Xs = filtfilt(b, a, double(X.')).';

    otherwise
        error('Unknown smoothing type: %s', cfg.type);
end
end

function w = getWin(cfg, hasFs, Fs, fld_ms, fld_samp)
w = [];
if hasFs && isfield(cfg, fld_ms) && ~isempty(cfg.(fld_ms))
    w = max(1, round(cfg.(fld_ms)/1000 * Fs));
elseif isfield(cfg, fld_samp) && ~isempty(cfg.(fld_samp))
    w = cfg.(fld_samp);
end
if isempty(w), w = 5; end   % sensible tiny fallback
end

function sig = getSigma(cfg, hasFs, Fs)
if hasFs && isfield(cfg,'sigma_ms') && ~isempty(cfg.sigma_ms)
    sig = max(1, cfg.sigma_ms/1000 * Fs);
elseif isfield(cfg,'sigma_samp') && ~isempty(cfg.sigma_samp)
    sig = cfg.sigma_samp;
else
    sig = 3;  % default few samples
end
end


