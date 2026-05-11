function run_spike_classifier_withTESTdata6(varargin)
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
% trials = build_trials_from_mat(filesS, filesN);

% save('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/Squiggles/func/DD_102025.mat', 'trials')
load('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/Squiggles/func/DD_102025.mat')

% Use these in the rest of your pipeline
data   = trials.data;        % {N×1} cell, each C×T
labels = trials.labels;      % categorical {'NoSpike','Spike'}
groups = trials.groups;      % {N×1} group IDs (leakage-safe split key)
refLbl = trials.refLbl;      % channel labels (row order of data epochs)

fprintf('Built trials: N=%d  |  C=%d  |  Spike=%d  |  NoSpike=%d\n', ...
    numel(data), numel(refLbl), sum(labels=='Spike'), sum(labels=='NoSpike'));

%% ------------------------ GROUPED TRAIN / VAL SPLIT ---------------------
[trainMask, valMask] = group_holdout_split(groups, opts.valFrac, opts.seed);
XTrain = data(trainMask);   YTrain = labels(trainMask);
XVal   = data(valMask);     YVal   = labels(valMask);

% Keep class order stable across Train/Val
cats   = categories(YTrain);
YTrain = reordercats(YTrain,cats);
YVal   = reordercats(YVal,cats);

assert(~isempty(XTrain), 'Training set is empty after group split. Reduce valFrac or check group IDs.');

%% ------------------------ Filtering (band-pass + notch) -------------------
fs = 500;
[Xtr, Xva] = filter_train_val(XTrain, XVal, fs, [4 30], 4, [], [], 30);

%% Per-epoch scaling
Xtr = maggrad_median_scale_per_epoch(Xtr, refLbl);
Xva = maggrad_median_scale_per_epoch(Xva, refLbl);

%% ---------- Augment TRAIN ONLY (end-to-end; no PCA needed) ----------
% % Build channel info (so spatial ops respect mag vs grad)
% chanInfo = struct();
% isMag  = endsWith(refLbl,'1');
% isGrad = endsWith(refLbl,'2') | endsWith(refLbl,'3');
% chanInfo.idxMag  = find(isMag);
% chanInfo.idxGrad = find(isGrad);
% 
% % Base augmentor you already have (crop + channel-shift + sign-flip + freq-mix)
% aug = default_meg_aug();          % your function
% aug.p_crop    = 0.85;             % more cropping variety
% aug.crop_frac_min = 0.6;          % keep at least 60% of the window
% aug.crop_frac_max = 1.0;
% aug.p_shift   = 0.5;  aug.max_shift = 8;    % gentle sensor-order shift (within mag/grad)
% aug.p_flip    = 0.5;                        % polarity flip
% aug.p_mix     = 0.4;  aug.mix_alpha_min = 0.15; aug.mix_alpha_max = 0.35;
% 
% % Negative pool for frequency-mix (NoSpike only)
% negPool = Xtr(YTrain=='NoSpike');
% 
% % Compose augmentor
% A = make_meg_augmentor(aug, negPool, chanInfo);
% 
% % Optional: add a few extra light ops (time-jitter, noise, channel-drop)
% jit_ms   = 8;    % ±8 ms shift
% jit_samp = max(1, round(jit_ms*fs/1000));
% p_noise  = 0.3;  % 30% add small Gaussian noise
% p_drop   = 0.15; % 15% random channel dropout (mask a few sensors)
% 
% Xtr_aug = cell(size(Xtr));
% for i = 1:numel(Xtr)
%     x = Xtr{i};
% 
%     % 0) (optional) small jitter around center (dont use circular if clips arent periodic)
%     if jit_samp>0
%         k = randi([-jit_samp, jit_samp],1,1);
%         s = max(1, 1+k); e = min(size(x,2), size(x,2)+k);
%         if k>=0, x = [x(:,1+k:end), x(:,end*ones(1,k))];
%         else,    x = [x(:,1*ones(1,-k)), x(:,1:end+k)];
%         end
%     end
% 
%     % 1) your compound MEG augmentor
%     x = A(x);
% 
%     % 2) light additive noise
%     if rand < p_noise
%         sigma = 0.02;                           % noise level after z-score
%         x = x + sigma*randn(size(x),'like',x);
%     end
% 
%     % 3) random channel dropout (simulate bad sensor / spatial robustness)
%     if rand < p_drop
%         nDrop = randi([1, max(1, round(0.03*size(x,1)))]);   % drop ~3% sensors
%         idxD  = randperm(size(x,1), nDrop);
%         x(idxD,:) = 0;
%     end
% 
%     Xtr_aug{i} = x;
% end
% 
% % Replace TRAIN with augmented set (you can also concatenate to enlarge TRAIN)
% Xtr = Xtr_aug;


%% Fit normalization on TRAIN only; apply to both
switch lower(opts.norm)
  case 'zscore'
    [mu,sd] = fit_channel_stats(Xtr);
    stats = struct('type','zscore','mu',mu(:),'sd',sd(:));
  case 'robust'
    [med,madv] = fit_channel_stats_robust(Xtr);
    stats = struct('type','robust','med',med(:),'madv',madv(:));
  otherwise, error('norm must be zscore or robust');
end
Xtr2 = normalize_cells(Xtr, stats, opts.norm);
Xva2 = normalize_cells(Xva, stats, opts.norm);

%% Fix a single length for both sets
Ttrain = cellfun(@(x) size(x,2), Xtr2);
Tval   = cellfun(@(x) size(x,2), Xva2);
Tfix   = max(1, round(prctile([Ttrain(:); Tval(:)], 95)));
Xtr2   = pad_sequences_to_length(Xtr2, Tfix);
Xva2   = pad_sequences_to_length(Xva2, Tfix);
padArgs = {"SequenceLength", Tfix, "SequencePaddingDirection","right"};

%% Class weights
cats = categories(YTrain);
YTrain = reordercats(YTrain,cats);  YVal = reordercats(YVal,cats);
cnt = countcats(YTrain); w = sum(cnt)./max(cnt,1); w = w/mean(w);
clsLayer = classificationLayer('Name','cls','Classes',cats,'ClassWeights',w);

%%
% Inputs you need
X = Xtr2{1442};                 % [C x T] epoch
fs = 500;                    % your sampling rate
ref = refLbl;                % channel labels matching rows of X

% pick a time window around peak GFP (change width as you like)
g = std(X,[],1); [~,ip] = max(g);
w_ms = 15; W = round(w_ms*fs/1000);             % e.g., 80 ms window
i1 = max(1, ip - floor(W/2)); i2 = min(size(X,2), ip + floor(W/2));
vals = mean(X(:, i1:i2), 2);                    % C×1 topography for the window

% build a minimal timelock struct
tl = [];
tl.label  = ref(:);
tl.dimord = 'chan_time';
tl.time   = 0;
tl.avg    = vals(:)';                           % 1×C (FieldTrip expects chan×time)

% choose a layout: mags vs planar grads
% If your epoch contains both: plot them separately.
cfg = [];
cfg.parameter = 'avg';
cfg.xlim      = [0 0];                          % single "time"
cfg.marker    = 'off'; cfg.comment = 'no';
% cfg.zlim = 'maxabs';  % or cfg.zlim = [-m m];
cfg.colorbar  = 'yes'; cfg.zlim = 'maxabs';

% Magnetometers
cfg.layout = 'neuromag306mag.lay';
ft_topoplotER(cfg, tl); title(sprintf('MEG topomap @ peak GFP (±%d ms)', w_ms));

%%
X = Xtr2{142};                 % [C x T] epoch
X = Xtr2{1142};                 % [C x T] epoch
X = Xtr2{1242};                 % [C x T] epoch
figure, plot(mean(mean(X),1))
g = std(X,[],1); [~,ip] = max(mean(mean(X),1)); ip = ip - 15;
w_ms = 10; W = round(w_ms*fs/1000);             % e.g., 80 ms window
i1 = max(1, ip - floor(W/2)); i2 = min(size(X,2), ip + floor(W/2));
vals = mean(X(:, i1:i2), 2);                    % C×1 topography for the window

% build a minimal timelock struct
tl = [];
tl.label  = ref(:);
tl.dimord = 'chan_time';
tl.time   = 0;
tl.avg    = vals(:)';                           % 1×C (FieldTrip expects chan×time)

% choose a layout: mags vs planar grads
% If your epoch contains both: plot them separately.
cfg = [];
cfg.parameter = 'avg';
cfg.xlim      = [0 0];                          % single "time"
cfg.marker    = 'off'; cfg.comment = 'no';
% cfg.zlim = 'maxabs';  % or cfg.zlim = [-m m];
cfg.colorbar  = 'yes'; cfg.zlim = 'maxabs';

% Magnetometers
cfg.layout = 'neuromag306mag.lay';
ft_topoplotER(cfg, tl); title(sprintf('MEG topomap @ peak GFP (±%d ms)', w_ms));

%%
% Xtr2 : {Ntrain×1} epochs (C×T), already filtered + mag/grad median scaled + normalized
% YTrain: categorical labels {'NoSpike','Spike'}
% fs    : 500
% refLbl: channel names (row order of Xtr2)
% raw_data.grad: sensor struct from your FT preprocessing (or sens)

[Topo, counts] = group_topo_mean(Xtr2, YTrain, fs, refLbl, 150, true);
Topo.counts = counts;  % optional, for titles

plot_group_topos_ft(Topo, raw_data.grad);


%%
% % Sanity: report sizes
% sz = cellfun(@(x) [size(x,1) size(x,2)], Xtr2, 'uni', false);
% rows = cellfun(@(s) s(1), sz); cols = cellfun(@(s) s(2), sz);
% fprintf('Rows unique: %s\n', mat2str(unique(rows)));
% fprintf('Cols range : [%d .. %d]\n', min(cols), max(cols));
% 
% % Re-pad with fixed Tfix using the corrected helper
% Xtr2 = pad_sequences_to_length(Xtr2, Tfix);
% Xva2 = pad_sequences_to_length(Xva2, Tfix);
% 
% % Remove any bad/empty epochs just in case
% C = size(Xtr2{1},1);
% goodTr = cellfun(@(x) isnumeric(x) && ismatrix(x) && size(x,1)==C && size(x,2)==Tfix && all(isfinite(x(:))), Xtr2);
% goodVa = cellfun(@(x) isnumeric(x) && ismatrix(x) && size(x,1)==C && size(x,2)==Tfix && all(isfinite(x(:))), Xva2);
% 
% if ~all(goodTr)
%     warning('Dropping %d malformed TRAIN epochs.', sum(~goodTr));
%     Xtr2 = Xtr2(goodTr);  YTrain = YTrain(goodTr);
% end
% if ~all(goodVa)
%     warning('Dropping %d malformed VAL epochs.', sum(~goodVa));
%     Xva2 = Xva2(goodVa);  YVal   = YVal(goodVa);
% end
% 
% % Ensure single precision
% Xtr2 = cellfun(@single, Xtr2, 'uni', false);
% Xva2 = cellfun(@single, Xva2, 'uni', false);

%% Sizes & environment
inputSize = size(Xtr2{1},1);   % number of sensors (C)
execEnv   = 'auto'; if opts.useGPU && gpuDeviceCount>0, execEnv='gpu'; end

%% Compact 1D CNN (temporal)
% after you compute Tfix and before training:
layers1d = [
    sequenceInputLayer(inputSize, ...
        'Normalization','none', ...
        'MinLength', Tfix, ...          % <= tell it your padded length
        'Name','in')

    convolution1dLayer(7, 32, 'Padding','same')
    reluLayer
    convolution1dLayer(5, 32, 'Padding','same')
    reluLayer
    maxPooling1dLayer(2,'Stride',2)     % requires MinLength >= 2
    convolution1dLayer(3, 64, 'Padding','same')
    reluLayer
    globalAveragePooling1dLayer
    dropoutLayer(0.2)
    fullyConnectedLayer(numel(cats))
    softmaxLayer
    classificationLayer('Name','cls','Classes',cats,'ClassWeights',w)
];

opts1d = trainingOptions('adam', ...
  'InitialLearnRate',3e-4, 'MaxEpochs',40, 'MiniBatchSize',256, ...
  padArgs{:}, ...
  'Shuffle','every-epoch', ...
  'ValidationData',{Xva2,YVal}, ...
  'ValidationFrequency',max(10,ceil(numel(Xtr2)/256)), ...
  'ValidationPatience',8, ...
  'ExecutionEnvironment',execEnv, ...
  'Plots','training-progress','Verbose',true);

net1d = trainNetwork(Xtr2, YTrain, layers1d, opts1d);

%% Pack cells to 4-D images [C x Tfix x 1 x N]
[ImTr, YimTr] = imgs4d(Xtr2, YTrain);   % helper below
[ImVa, YimVa] = imgs4d(Xva2, YVal);

layers2d = [
    imageInputLayer([size(ImTr,1) size(ImTr,2) 1], 'Normalization','none', 'Name','in')

    % temporal-focused convs
    convolution2dLayer([1 7], 32, 'Padding','same', 'Name','tconv7')
    reluLayer('Name','trelu7')
    convolution2dLayer([1 5], 32, 'Padding','same', 'Name','tconv5')
    reluLayer('Name','trelu5')

    % spatial mixing across sensors (rows)
    convolution2dLayer([5 1], 32, 'Padding','same', 'Name','sconv5')
    reluLayer('Name','srelu5')

    % joint spatio-temporal
    convolution2dLayer([3 3], 64, 'Padding','same', 'Name','st33')
    reluLayer('Name','strelu33')
    maxPooling2dLayer([2 2], 'Stride',[2 2], 'Name','pool1')

    convolution2dLayer([3 3], 128, 'Padding','same', 'Name','c3')
    reluLayer('Name','r3')
    globalAveragePooling2dLayer('Name','gap')

    dropoutLayer(0.3, 'Name','drop')
    fullyConnectedLayer(numel(cats), 'Name','fc')
    softmaxLayer('Name','sm')
    classificationLayer('Name','cls', 'Classes',cats, 'ClassWeights',w)
];

opts2d = trainingOptions('adam', ...
    'InitialLearnRate',3e-4, 'MaxEpochs',40, 'MiniBatchSize',128, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{ImVa,YimVa}, ...
    'ValidationFrequency',max(10,ceil(size(ImTr,4)/128)), ...
    'ValidationPatience',8, ...
    'ExecutionEnvironment',execEnv, ...
    'Plots','training-progress','Verbose',true);

net2d = trainNetwork(ImTr, YimTr, layers2d, opts2d);

%%
[cats] = categories(YVal);
[alpha, BAval, cmVal] = pick_fusion_alpha(net1d, net2d, Xva2, ImVa, YVal, cats);

%%
% Build fusion features on VAL
[~, s1v] = classify(net1d, Xva2);  p1v = s1v(:, strcmp(cats,'Spike'));
[~, s2v] = classify(net2d, ImVa);  p2v = s2v(:, strcmp(cats,'Spike'));
Fv = [p1v p2v];

% Train logistic meta on VAL (or CV folds if you want)
meta = fitclinear(Fv, YVal, 'Learner','logistic', 'Solver','lbfgs', 'ClassNames',cats);

% Predict on TEST (or on VAL for check)
% Build Ftest = [p1test p2test] the same way, then:
Yp_meta = predict(meta, Fv);
cm  = confusionmat(YVal, Yp_meta, 'Order', cats);
rec = diag(cm)./max(1,sum(cm,2)); BA = mean(rec);
fprintf('Stacking BA=%.3f\n', BA);


%%
% Ensure identical categorical types & order
cats   = categories(YVal);                            % reference order
YvalU  = categorical(YVal, cats, 'Ordinal', false);   % force order

% Coerce Yp_meta to the same categorical
if isa(Yp_meta,'categorical')
    YpU = reordercats(Yp_meta, cats);
elseif isstring(Yp_meta) || iscellstr(Yp_meta)
    YpU = categorical(Yp_meta, cats, 'Ordinal', false);
elseif isnumeric(Yp_meta)
    % if numeric 1..K, map to cats
    YpU = categorical(Yp_meta, 1:numel(cats), cats, 'Ordinal', false);
else
    error('Unexpected type for Yp_meta: %s', class(Yp_meta));
end

% If the predictor produced any unknown labels, merge them to a known class (e.g., majority)
extra = setdiff(categories(YpU), cats);
if ~isempty(extra)
    YpU = mergecats(YpU, extra, cats{1});
    YpU = reordercats(YpU, cats);
end

% Now safe:
cm  = confusionmat(YvalU, YpU, 'Order', cats);
acc = sum(diag(cm))/sum(cm,'all');
rec = diag(cm)./max(1, sum(cm,2));
BA  = mean(rec);

fprintf('Acc=%.3f  BA=%.3f\n', acc, BA);





































%%
cfg = []; cfg.layout = 'neuromag306mag.lay'; lay = ft_prepare_layout(cfg);


filename = '/data/MEG/Research/SpikeDetection/tsss_ecgClean_Data/alioto_courtney_Run03_spont_supine_raw_t_sss_ecgClean_raw_DS.fif';
cfg = []; cfg.dataset = filename;  %
cfg.channel = {'meg'};
raw_data = ft_preprocessing(cfg);

chlabels=raw_data.label;
isMag  = endsWith(chlabels,'1');
isGrad  = endsWith(chlabels,'2') | endsWith(chlabels,'3');

%% Smoothing
smoothing = 0;
if smoothing == 1
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
end

%% ----------------------- CHANNEL-WISE NORMALIZATION ---------------------
switch lower(opts.norm)
    case 'zscore'
        [mu, sd] = fit_channel_stats(Xtr);
        mu = ensure_column(mu); sd = ensure_column(sd);
        zscoref  = @(x) bsxfun(@rdivide, bsxfun(@minus,x,mu), max(sd,1e-12));
        normStats = struct('type','zscore','mu',mu,'sd',sd);
    case 'robust'
        [med, madv] = fit_channel_stats_robust(Xtr);
        med = ensure_column(med); madv = ensure_column(madv);
        zscoref  = @(x) bsxfun(@rdivide, bsxfun(@minus,x,med), max(madv,1e-12));
        normStats = struct('type','robust','med',med,'madv',madv);
        mu = []; sd = [];
    otherwise
        error('Unknown opts.norm: %s', opts.norm);
end
XTrain2 = cellfun(zscoref, Xtr,'UniformOutput',false);
XVal2   = cellfun(zscoref, XVal  ,'UniformOutput',false);

%% ----- Fix a single SequenceLength for both sets (even when PCA is off)
% opts.pca = 0;
% if opts.pca == 1
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
% end

execEnv = 'auto';
if opts.useGPU && gpuDeviceCount>0, execEnv = 'gpu'; end
valFreq = max(5, ceil(numel(XTrain2)/128));

switch lower(opts.arch)
    % ---------- 1) ResNet1D (sequence; uses your helper) ----------
    case 'resnet1d'
        lgraph = resnet1d_layers(inputSize, numel(cats));
        lgraph = replaceLayer(lgraph, "cls", clsLayer);  % class weights

        trainOpts = trainingOptions('adam', ...
            'InitialLearnRate',3e-4, 'MaxEpochs',60, 'MiniBatchSize',128, ...
            padArgs{:}, ...
            'Shuffle','every-epoch', ...
            'ValidationData',{XVal2,YVal}, ...
            'ValidationFrequency',valFreq, 'ValidationPatience',10, ...
            'ExecutionEnvironment',execEnv, ...
            'Verbose',opts.verbose, 'Plots','training-progress');

        net = trainNetwork(XTrain2, YTrain, lgraph, trainOpts);
        mdlName = 'resnet1d_spike_classifier.mat';

    % ---------- 2) Mamba-proxy (dilated Conv1D + GLU + GRU) ----------
    case 'mamba_proxy'
        lgraph = mamba_proxy_layers(inputSize, numel(cats), clsLayer);  % from earlier messages
        trainOpts = trainingOptions('sgdm', ...
            'InitialLearnRate',5e-3, 'Momentum',0.9, 'MaxEpochs',20, 'MiniBatchSize',256, ...
            padArgs{:}, ...
            'GradientThresholdMethod','l2norm','GradientThreshold',1, ...
            'Shuffle','every-epoch', ...
            'ValidationData',{XVal2,YVal}, ...
            'ValidationFrequency',max(20,ceil(numel(XTrain2)/256)), 'ValidationPatience',5, ...
            'ExecutionEnvironment',execEnv, 'Plots','training-progress', 'Verbose',opts.verbose);

        net = trainNetwork(XTrain2, YTrain, lgraph, trainOpts);
        mdlName = 'mamba_proxy_spike_classifier.mat';

    % ---------- 3) LV-CadeNet(paper-style) (image input; requires LV stacks) ----------
    case 'lv_cadenet_paper'
        % Build LV stacks ONCE (raw+6 long-view maps), pack to [C×L×7×N]
        fs = 200;                                % set your true post-downsample fs
        XTrLV = cellfun(@(x) build_longview_features(single(x), fs), XTrain2, 'uni', false);
        XVaLV = cellfun(@(x) build_longview_features(single(x), fs), XVal2,   'uni', false);
        L = max(1, round(prctile(cellfun(@(x) size(x,2), XTrLV), 85)));
        [trainImgs, trainLabs] = cell2imgs(XTrLV, YTrain, L);  % [C×L×7×N]
        [valImgs,   valLabs]   = cell2imgs(XVaLV, YVal,   L);

        C = size(trainImgs,1); D = size(trainImgs,3);
        lgraph = lv_cadenet_paper_layers(C, L, D, numel(cats), clsLayer);

        trainOpts = trainingOptions('adam', ...
            'InitialLearnRate',3e-4, 'MaxEpochs',40, 'MiniBatchSize',128, ...
            'Shuffle','every-epoch', ...
            'ValidationData',{valImgs,valLabs}, ...
            'ValidationFrequency',max(10,ceil(size(trainImgs,4)/128)), 'ValidationPatience',8, ...
            'L2Regularization',1e-4, ...
            'ExecutionEnvironment',execEnv, 'Plots','training-progress', 'Verbose',opts.verbose);

        net = trainNetwork(trainImgs, trainLabs, lgraph, trainOpts);
        mdlName = 'lv_cadenet_paper.mat';

    % ---------- 4) Tiny LSTM (baseline) ----------
    case 'lstm'
        layers = [
          sequenceInputLayer(inputSize,'Normalization','none')
          bilstmLayer(64,'OutputMode','last')
          dropoutLayer(0.2)
          fullyConnectedLayer(numel(cats))
          softmaxLayer
          classificationLayer('Name','cls','Classes',cats,'ClassWeights',w)];
        trainOpts = trainingOptions('adam', ...
            'InitialLearnRate',1e-3, 'MaxEpochs',40, 'MiniBatchSize',128, ...
            padArgs{:}, ...
            'Shuffle','every-epoch', ...
            'ValidationData',{XVal2,YVal}, ...
            'ValidationFrequency',valFreq, 'ValidationPatience',10, ...
            'ExecutionEnvironment',execEnv, 'Plots','training-progress', 'Verbose',opts.verbose);
        net = trainNetwork(XTrain2, YTrain, layers, trainOpts);
        mdlName = 'lstm_spike_classifier.mat';

    otherwise
        error('Unknown arch: %s. Use ''resnet1d'', ''mamba_proxy'', ''lv_cadenet_paper'', or ''lstm''.',opts.arch);
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


