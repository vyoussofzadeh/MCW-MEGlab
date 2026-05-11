function run_spike_classifier_withTESTdata9(varargin)
%% Spike-Detection MEG Pipeline  Spike vs NoSpike (clean + dual-branch)
% Author: MCW MEG Lab  V. Youssofzadeh <vyoussofzadeh@mcw.edu>

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

% Data handling
opts.capPerRun   = 0;                  % 0=off; else cap epochs per class per run (e.g., 200)

% Preprocessing
opts.fs          = 500;                % Hz
opts.bpHz        = [4 40];             % band-pass
opts.bpOrder     = 4;                  % IIR order
opts.notchHz     = [];                 % [] to disable (e.g., 50 or 60 to enable)
opts.notchQ      = 35;                 % quality factor (IIR notch), used if notchHz set

% Normalization
opts.norm        = 'robust';           % 'robust' | 'zscore'

% Temporal window and padding
opts.tfixPct     = 80;                 % percentile for fixed length padding (shorter = crisper spikes)

% Augmentation (TRAIN only)
opts.useAug      = true;               % turn on moderate, MEG-aware augs
opts.jitter_ms   = 8;                  % ± jitter along time
opts.p_noise     = 0.25;               % additive noise probability
opts.noiseSigma  = 0.02;               % z-units after norm
opts.p_drop      = 0.10;               % random channel dropout prob
opts.dropFrac    = 0.03;               % fraction of channels to zero out

% Temporal branch PCA
opts.usePCA      = true;               % PCA across channels for temporal models
opts.varKeep     = 98;                 % % variance to keep for PCA (9899 recommended)

% Models
opts.trainLSTM   = true;
opts.trainCNN1D  = true;
opts.trainRes1D  = true;
opts.trainCNN2D  = true;

% Fusion
opts.doFusion    = true;

% RNG
opts.seed        = 0;
rng(opts.seed,'twister');

% Allow name-value overrides
opts = parse_inputs(opts,varargin{:});

%% ----------------------------- ENV SETUP --------------------------------
setup_environment(opts);
if ~exist(opts.modelDir,'dir'); mkdir(opts.modelDir); end
execEnv = 'auto';
if opts.useGPU
    try, g=gpuDevice; if ~isempty(g)&&g.SupportsDouble, execEnv='gpu'; end, catch, execEnv='auto'; end
end

%% ------------------------- LOAD / BUILD DATA ----------------------------
% If you prefer to build from files, uncomment the two lines below:
filesS = list_mat_files(opts.spikeDir);
filesN = list_mat_files(opts.noSpikeDir);
% trials = build_trials_from_mat(filesS, filesN);

% Using your cached trials struct:
load('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/Squiggles/func/DD_102025.mat');

data   = trials.data;        % {N×1}, each C×T
labels = trials.labels;      % categorical {'NoSpike','Spike'}
groups = trials.groups;      % {N×1} run/subject IDs
refLbl = trials.refLbl;      % C×1 channel names

times = linspace(-200,200,size(data{1},2));

fprintf('Trials: N=%d | C=%d | Spike=%d | NoSpike=%d\n', ...
    numel(data), numel(refLbl), sum(labels=='Spike'), sum(labels=='NoSpike'));

%%
% 1) counts
fprintf('N=%d | Spike=%d | NoSpike=%d | C=%d\n', numel(trials.data), sum(trials.labels=='Spike'), sum(trials.labels=='NoSpike'), numel(trials.refLbl));

% 2) epoch size consistency
sz = cellfun(@size, trials.data, 'uni', false);
rows = cellfun(@(s) s(1), sz); cols = cellfun(@(s) s(2), sz);
disp(['Unique C (must be 1 value): ' mat2str(unique(rows))]);
% figure; histogram(cols, 40); title('T per epoch'); grid on

% 3) quick GFP plots
% I = randperm(numel(trials.data), min(6,numel(trials.data)));
% figure('Color','w'); tiledlayout(numel(I),1,'TileSpacing','compact');
% for k=1:numel(I), g=std(trials.data{I(k)},[],1); nexttile; plot(g,'k'); grid on; title(sprintf('epoch %d', I(k))); end

%% ------------------------ CAP BY RUN (optional) -------------------------
if opts.capPerRun>0
    [data, labels, groups] = cap_by_group(data, labels, groups, opts.capPerRun);
    fprintf('After capPerRun=%d --> N=%d\n', opts.capPerRun, numel(data));
end

%% ------------------------ GROUPED TRAIN / VAL SPLIT ---------------------
[trainMask, valMask] = group_holdout_split(groups, opts.valFrac, opts.seed);
XTrain = data(trainMask);   YTrain = labels(trainMask);
XVal   = data(valMask);     YVal   = labels(valMask);

cats   = categories(YTrain);
YTrain = reordercats(YTrain,cats);
YVal   = reordercats(YVal,cats);
assert(~isempty(XTrain), 'Empty TRAIN after group split.');

fprintf('TRAIN: Spike=%d, NoSpike=%d | VAL: Spike=%d, NoSpike=%d | Groups train=%d, val=%d\n', ...
    sum(YTrain=='Spike'), sum(YTrain=='NoSpike'), sum(YVal=='Spike'), sum(YVal=='NoSpike'), ...
    numel(unique(groups(trainMask))), numel(unique(groups(valMask))));

%% --------------------------- FILTERING ----------------------------------
[Xtr, Xva] = filter_train_val(XTrain, XVal, opts.fs, opts.bpHz, opts.bpOrder, opts.notchHz, opts.notchQ);

%% -------------- MAG/GRAD MEDIAN SCALING (per epoch) --------------------
Xtr = maggrad_median_scale_per_epoch(Xtr, refLbl);
Xva = maggrad_median_scale_per_epoch(Xva, refLbl);

%% -------------------------- AUGMENT (TRAIN) -----------------------------
if opts.useAug
    Xtr = augment_meg_train(Xtr, YTrain, refLbl, opts);
end

%% ----------------------- NORMALIZATION (fit on TRAIN) -------------------
switch lower(opts.norm)
    case 'zscore'
        [mu,sd] = fit_channel_stats(Xtr);
        stats = struct('type','zscore','mu',mu(:),'sd',sd(:));
    case 'robust'
        [med,madv] = fit_channel_stats_robust(Xtr);
        stats = struct('type','robust','med',med(:),'madv',madv(:));
    otherwise, error('opts.norm must be zscore or robust');
end
Xtr2 = normalize_cells(Xtr, stats);
Xva2 = normalize_cells(Xva, stats);

%% ----------------------- FIXED LENGTH (shorter) -------------------------
Ttrain = cellfun(@(x) size(x,2), Xtr2);
Tval   = cellfun(@(x) size(x,2), Xva2);
Tfix   = max(1, round(prctile([Ttrain(:); Tval(:)], opts.tfixPct)));

Xtr2 = warp_sequences_to_length_interp(Xtr2, Tfix, 'pchip');  % Option B
Xva2 = warp_sequences_to_length_interp(Xva2, Tfix, 'pchip');
padArgs = {"SequenceLength", Tfix, "SequencePaddingDirection","right"};

fprintf('Tfix=%d (percentile %d)\n', Tfix, opts.tfixPct);

%% --------------------------- CLASS WEIGHTS ------------------------------
cnt = countcats(YTrain); w = sum(cnt)./max(cnt,1); w = w/mean(w);
clsLayer = classificationLayer('Name','cls','Classes',cats,'ClassWeights',w);

%% ------------------------ BRANCH A: TEMPORAL ----------------------------
% PCA for temporal models ONLY (keep 2D branch in sensor space)

if opts.usePCA
    [W_pca, mC_pca, expl, K] = pca_fit_channels_cell(Xtr2, opts.varKeep);
    Xtr_tem = pca_apply_channels_cell(Xtr2, W_pca, mC_pca, K);
    Xva_tem = pca_apply_channels_cell(Xva2, W_pca, mC_pca, K);
    inputSize_tem = K;
    fprintf('PCA (temporal): K=%d (%.1f%% var)\n', K, sum(expl(1:K)));
else
    Xtr_tem = Xtr2; Xva_tem = Xva2; inputSize_tem = size(Xtr2{1},1);
end

%%
%% ==== RMS representation for 2D CNN branch ====

win_sec  = 0.020;                        % 20 ms RMS window (tune this)
win_samp = max(1, round(win_sec * opts.fs));

% TRAIN: compute RMS over time for each epoch
Xtr_rms = cell(size(Xtr2));
for i = 1:numel(Xtr2)
    d = Xtr2{i};                         % [C x Tfix]
    Xtr_rms{i} = sqrt(movmean(d.^2, win_samp, 2));   % [C x Tfix]
end

% VAL: same
Xva_rms = cell(size(Xva2));
for i = 1:numel(Xva2)
    d = Xva2{i};
    Xva_rms{i} = sqrt(movmean(d.^2, win_samp, 2));   % [C x Tfix]
end

% Convert to 4D images [C x T x 1 x N]
[ImTr_rms, YTr_rms] = imgs4d(Xtr_rms, YTrain);   % C x T x 1 x N
[ImVa_rms, YVa_rms] = imgs4d(Xva_rms, YVal);     % C x T x 1 x N

C = size(ImTr_rms,1);
T = size(ImTr_rms,2);
D = size(ImTr_rms,3);   % should be 1

%%
% Instead of Xtr_tem = Xtr2; you can do:
Xtr_tem = Xtr_rms;
Xva_tem = Xva_rms;
inputSize_tem = size(Xtr_rms{1},1);   % still C

%% ----------------------------- MODELS -----------------------------------
models = struct();

% --- LSTM ---
if opts.trainLSTM

    layers_lstm = [
        sequenceInputLayer(inputSize_tem,'Normalization','none','MinLength',Tfix,'Name','in')
        bilstmLayer(64,'OutputMode','last','Name','bilstm')
        dropoutLayer(0.2,'Name','drop')
        fullyConnectedLayer(numel(cats),'Name','fc')
        softmaxLayer('Name','sm')
        classificationLayer('Name','cls','Classes',cats,'ClassWeights',w)];
    opts_lstm = trainingOptions('adam', ...
        'InitialLearnRate',3e-4,'MaxEpochs',40,'MiniBatchSize',256, ...
        padArgs{:}, 'Shuffle','every-epoch', ...
        'ValidationData',{Xva_tem,YVal}, ...
        'ValidationFrequency',max(10,ceil(numel(Xtr_tem)/256)), ...
        'ValidationPatience',8, 'ExecutionEnvironment',execEnv, ...
        'Plots','training-progress','Verbose',true);
    models.net_lstm = trainNetwork(Xtr_tem, YTrain, layers_lstm, opts_lstm);
    [Yhat, s] = classify(models.net_lstm, Xva_tem);
    [acc, BA, cm] = metrics(YVal, Yhat); fprintf('[LSTM] Acc=%.3f BA=%.3f\n',acc,BA);
    models.lstm_val = struct('acc',acc,'BA',BA,'cm',cm,'scores',s);

end

% --- 1D CNN (compact) ---
if opts.trainCNN1D

    % Training options expect fixed sequence length:
    padArgs = {'SequenceLength', Tsmall, 'SequencePaddingDirection','right'};

    layers1d = [
        sequenceInputLayer(inputSize_tem,'Normalization','none','MinLength',Tsmall,'Name','in')
        convolution1dLayer(5,64,'Padding','same')
        batchNormalizationLayer
        reluLayer
        convolution1dLayer(3,32,'Padding','same')
        batchNormalizationLayer
        reluLayer
        convolution1dLayer(3,64,'Padding','same')
        batchNormalizationLayer
        reluLayer
        globalAveragePooling1dLayer
        dropoutLayer(0.3)
        fullyConnectedLayer(128)
        reluLayer
        dropoutLayer(0.3)
        fullyConnectedLayer(numel(cats))
        softmaxLayer
        classificationLayer('Name','cls','Classes',cats,'ClassWeights',w)
        ];

    opts1d = trainingOptions('adam', ...
        'InitialLearnRate',3e-4,'MaxEpochs',60,'MiniBatchSize',256, ...
        'L2Regularization',1e-4,'GradientThreshold',1, ...
        padArgs{:}, 'Shuffle','every-epoch', ...
        'ValidationData',{Xva_tem,YVal}, ...
        'ValidationFrequency',max(10,ceil(numel(Xtr_tem)/256)), ...
        'ValidationPatience',8,'ExecutionEnvironment',execEnv, ...
        'Plots','training-progress','Verbose',true);

    net1d = trainNetwork(Xtr_tem, YTrain, layers1d, opts1d);
    [Yhat1d, s1d] = classify(net1d, Xva_tem);

end

% --- ResNet-1D (dilated) ---
if opts.trainRes1D

    lgraph = resnet1d_layers(inputSize_tem, numel(cats));
    lgraph = replaceLayer(lgraph, "cls", clsLayer);
    mb = 128; lr = 3e-4; maxE = 120;
    valFreq = max(5, ceil(numel(Xtr_tem)/mb));
    optsR = trainingOptions('adam', ...
        'InitialLearnRate',lr,'LearnRateSchedule','piecewise',...
        'LearnRateDropPeriod',25,'LearnRateDropFactor',0.5, ...
        'MaxEpochs',maxE,'MiniBatchSize',mb, 'L2Regularization',5e-4, ...
        'GradientThresholdMethod','l2norm','GradientThreshold',1, ...
        padArgs{:}, 'Shuffle','every-epoch', ...
        'ValidationData',{Xva_tem,YVal}, 'ValidationFrequency',valFreq, ...
        'ValidationPatience',10, 'ExecutionEnvironment',execEnv, ...
        'Plots','training-progress','Verbose',true,'CheckpointPath',opts.modelDir);
    models.net_res1d = trainNetwork(Xtr_tem, YTrain, lgraph, optsR);
    [Yhat, s] = classify(models.net_res1d, Xva_tem);
    [acc, BA, cm] = metrics(YVal, Yhat); fprintf('[ResNet1D] Acc=%.3f BA=%.3f\n',acc,BA);
    models.res1d_val = struct('acc',acc,'BA',BA,'cm',cm,'scores',s);

end

% --- 2D CNN (sensor space) ---
if opts.trainCNN2D

    C = size(ImTr_rms,1); 
    T = size(ImTr_rms,2); 
    D = size(ImTr_rms,3);   % 1

    % Per-channel normalization (fit on TRAIN RMS only)
    Xall = reshape(ImTr_rms, C, []);          % collapse (T*D*batch)
    mu = mean(Xall,2); 
    sd = std(Xall,0,2);  sd(sd<1e-12)=1;

    for i=1:size(ImTr_rms,4)
        for d=1:D
            ImTr_rms(:,:,d,i) = (ImTr_rms(:,:,d,i) - mu) ./ sd;
        end
    end
    for i=1:size(ImVa_rms,4)
        for d=1:D
            ImVa_rms(:,:,d,i) = (ImVa_rms(:,:,d,i) - mu) ./ sd;
        end
    end

    % imageInputLayer expects [H W C] = [C T D]
    layers2d = [
        imageInputLayer([C T D],'Normalization','none','Name','in')
        convolution2dLayer([1 7],32,'Padding','same')
        batchNormalizationLayer
        reluLayer
        convolution2dLayer([1 5],32,'Padding','same')
        batchNormalizationLayer
        reluLayer
        convolution2dLayer([5 1],32,'Padding','same')
        batchNormalizationLayer
        reluLayer
        convolution2dLayer([3 3],64,'Padding','same')
        batchNormalizationLayer
        reluLayer
        maxPooling2dLayer([2 2],'Stride',[2 2])
        convolution2dLayer([3 3],128,'Padding','same')
        batchNormalizationLayer
        reluLayer
        globalAveragePooling2dLayer
        dropoutLayer(0.3)
        fullyConnectedLayer(numel(cats))
        softmaxLayer
        classificationLayer('Name','cls','Classes',cats,'ClassWeights',w)
        ];

    valFreq = max(10, ceil(size(ImTr_rms,4)/128));
    opts2d = trainingOptions('adam', ...
        'InitialLearnRate',3e-4,'MaxEpochs',60,'MiniBatchSize',128, ...
        'L2Regularization',1e-4, ...
        'LearnRateSchedule','piecewise','LearnRateDropFactor',0.5,'LearnRateDropPeriod',8, ...
        'Shuffle','every-epoch', ...
        'ValidationData',{ImVa_rms,YVa_rms}, ...
        'ValidationFrequency',valFreq, ...
        'ValidationPatience',8, ...
        'ExecutionEnvironment',execEnv, ...
        'Plots','training-progress','Verbose',true);

    models.net_2d_rms = trainNetwork(ImTr_rms, YTr_rms, layers2d, opts2d);

    % Evaluate
    [Yhat2d_rms, s2d_rms] = classify(models.net_2d_rms, ImVa_rms);
    cm2d = confusionmat(YVa_rms, Yhat2d_rms, 'Order', cats);
    acc2d = sum(diag(cm2d))/sum(cm2d,'all');
    rec2d = diag(cm2d)./max(1,sum(cm2d,2));
    BA2d  = mean(rec2d);
    fprintf('[CNN2D-RMS] Acc=%.3f BA=%.3f\n', acc2d, BA2d);
    models.cnn2d_rms_val = struct('acc',acc2d,'BA',BA2d,'cm',cm2d,'scores',s2d_rms,'ImVa',ImVa_rms);

end

% --- 2D CNN (sensor space) ---
if opts.trainCNN2D

    C = size(ImTr,1); Tfix = size(ImTr,2); D = size(ImTr,3);

    layers2d = [
        imageInputLayer([size(ImTr,1) size(ImTr,2) size(ImTr,3)], 'Normalization','none')
        convolution2dLayer([1 7], 32, 'Padding','same')
        reluLayer
        convolution2dLayer([1 5], 32, 'Padding','same')
        reluLayer
        convolution2dLayer([5 1], 32, 'Padding','same')
        reluLayer
        convolution2dLayer([3 3], 64, 'Padding','same')
        reluLayer
        maxPooling2dLayer([2 2],'Stride',[2 2])
        convolution2dLayer([3 3], 128, 'Padding','same')
        reluLayer
        globalAveragePooling2dLayer
        dropoutLayer(0.3)
        fullyConnectedLayer(numel(cats))
        softmaxLayer
        classificationLayer('Name','cls','Classes',cats,'ClassWeights',w)
        ];

    opts2d = trainingOptions('adam', ...
        'InitialLearnRate',3e-4, 'MaxEpochs',40, 'MiniBatchSize',128, ...
        'Shuffle','every-epoch', ...
        'ValidationData',{ImVa,YimVa}, ...
        'ValidationFrequency', max(10, ceil(size(ImTr,4)/128)), ...
        'ValidationPatience',8, ...
        'ExecutionEnvironment',execEnv, ...
        'Plots','training-progress','Verbose',true);

    net2d = trainNetwork(ImTr, YimTr, layers2d, opts2d);

end

%% ----------------------------- FUSION -----------------------------------
if opts.doFusion && isfield(models,'net_2d')
    % pick p(spike) from best temporal net if available
    if isfield(models,'net_res1d'), netA = models.net_res1d; XvaA = Xva_tem;
    elseif isfield(models,'net_1d'), netA = models.net_1d;  XvaA = Xva_tem;
    elseif isfield(models,'net_lstm'), netA = models.net_lstm; XvaA = Xva_tem;
    else, netA = []; end

    if ~isempty(netA)
        [alpha, BAval, cmVal] = pick_fusion_alpha(netA, models.net_2d, XvaA, models.cnn2d_val.ImVa, YVal, cats);
        fprintf('[Fusion] alpha=%.2f | BA=%.3f\n', alpha, BAval);
        models.fusion = struct('alpha',alpha,'BA',BAval,'cm',cmVal);
    end
end

end % ========================= END MAIN ===================================


%% ========================= HELPER FUNCTIONS =============================

function setup_environment(o)
addpath(o.ftRoot); ft_defaults;
addpath(genpath(fullfile(o.codeRoot,'FT_functions','functions_new')));
addpath(genpath(fullfile(o.codeRoot,'FT_functions','helper')));
end

function files = list_mat_files(root)
if exist('rdir','file')==2
    files = rdir(fullfile(root,'*.mat'));
else
    dd = dir(fullfile(root,'**','*.mat'));
    for i=1:numel(dd), dd(i).name = fullfile(dd(i).folder, dd(i).name); end
    files = dd;
end
end

function [trainMask,valMask] = group_holdout_split(groups, valFrac, seed)
if nargin<3, seed=0; end, rng(seed,'twister');
u = unique(groups,'stable'); G=numel(u);
if G==0, error('No groups'); end
if G==1, valMask=false(size(groups)); trainMask=true(size(groups)); return; end
nVal = max(1, min(round(valFrac*G), G-1));
idx = randperm(G); valGroups = u(idx(1:nVal));
valMask   = ismember(groups, valGroups);
trainMask = ~valMask;
end

function [Xtr_f, Xva_f] = filter_train_val(Xtr, Xva, fs, bpHz, bpOrder, notchHz, notchQ)
% Robust IIR zero-phase filtering (Butter); optional notch
[Xtr_f, Xva_f] = deal(Xtr, Xva);
[b,a] = butter(bpOrder, bpHz/(fs/2), 'bandpass');
useNotch = ~isempty(notchHz) && isnumeric(notchHz) && isfinite(notchHz) && notchHz>0 && notchHz<fs/2;
if useNotch
    wo = notchHz/(fs/2); bw = wo/notchQ;
    [bn,an] = iirnotch(wo, bw);
end
filtf = @(x) (filtfilt(b,a,double(x.')).');   % (C×T)
for i=1:numel(Xtr_f)
    xi = Xtr_f{i}; xi = filtf(xi);
    if useNotch, xi = (filtfilt(bn,an,double(xi.')).'); end
    Xtr_f{i} = xi;
end
for i=1:numel(Xva_f)
    xi = Xva_f{i}; xi = filtf(xi);
    if useNotch, xi = (filtfilt(bn,an,double(xi.')).'); end
    Xva_f{i} = xi;
end
end

function Xs = maggrad_median_scale_per_epoch(Xcells, refLbl)
isMag  = endsWith(refLbl,'1');
isGrad = endsWith(refLbl,'2') | endsWith(refLbl,'3');
Xs = cell(size(Xcells));
for i=1:numel(Xcells)
    X = Xcells{i};
    if any(isMag)
        mMag  = median(abs(X(isMag,:)),'all','omitnan') + eps;
        X(isMag,:)  = X(isMag,:) / mMag;
    end
    if any(isGrad)
        mGrad = median(abs(X(isGrad,:)),'all','omitnan') + eps;
        X(isGrad,:) = X(isGrad,:) / mGrad;
    end
    Xs{i} = X;
end
end

function Xaug = augment_meg_train(Xcells, Y, refLbl, o)
% Moderate, MEG-aware augmentation
isMag  = endsWith(refLbl,'1');
isGrad = endsWith(refLbl,'2') | endsWith(refLbl,'3');
negPool = Xcells(Y=='NoSpike');
Xaug = cell(size(Xcells));

for i=1:numel(Xcells)
    x = Xcells{i};
    % Jitter (± o.jitter_ms)
    if o.jitter_ms>0
        J = max(1,round(o.jitter_ms*o.fs/1000));
        k = randi([-J,J],1,1);
        if k>0, x = [x(:,1+k:end) x(:,end*ones(1,k))];
        elseif k<0, x = [x(:,1*ones(1,-k)) x(:,1:end+k)];
        end
    end
    % Small crop (0.7..1.0)
    frac = 0.7 + 0.3*rand;
    T = size(x,2); w = max(1,round(frac*T)); s = randi([1 max(1,T-w+1)]); x = x(:, s:s+w-1);

    % Polarity flip
    if rand < 0.5, x = -x; end

    % Frequency mix with a random negative
    if ~isempty(negPool) && rand < 0.3
        p = negPool{randi(numel(negPool))};
        Tm = min(size(x,2), size(p,2));
        X = fft(x(:,1:Tm),[],2); P = fft(p(:,1:Tm),[],2);
        alpha = 0.15 + 0.2*rand;
        mag = (1-alpha)*abs(X) + alpha*abs(P);
        x(:,1:Tm) = real(ifft(mag.*exp(1j*angle(X)),[],2));
    end

    % Light noise
    if rand < o.p_noise
        x = x + o.noiseSigma * randn(size(x),'like',x);
    end

    % Random channel dropout (few sensors)
    if rand < o.p_drop
        nDrop = max(1, round(o.dropFrac*size(x,1)));
        idxD = randperm(size(x,1), nDrop);
        x(idxD,:) = 0;
    end

    Xaug{i} = x;
end
end

function Xn = normalize_cells(Xcells, stats)
Xn = cell(size(Xcells));
switch lower(stats.type)
    case 'zscore'
        mu = stats.mu(:); sd = max(stats.sd(:),1e-12);
        for i=1:numel(Xcells)
            x = Xcells{i}; Xn{i} = (x - mu) ./ sd;
        end
    case 'robust'
        med = stats.med(:); madv = max(stats.madv(:),1e-12);
        for i=1:numel(Xcells)
            x = Xcells{i}; Xn{i} = (x - med) ./ madv;
        end
    otherwise, error('stats.type must be zscore or robust');
end
end

function Xp = pad_sequences_to_length(Xcells, Tfix)
Xp = cell(size(Xcells));
for i=1:numel(Xcells)
    x = Xcells{i};
    if size(x,2) >= Tfix
        Xp{i} = x(:,1:Tfix);
    else
        pad = repmat(x(:,end), 1, Tfix - size(x,2));
        Xp{i} = [x pad];
    end
end
end

function [Im4d, Y1] = imgs4d(Xcells, Y)
N = numel(Xcells); C = size(Xcells{1},1); T = size(Xcells{1},2);
Im4d = zeros(C, T, 1, N, 'single');
for i=1:N
    Xi = Xcells{i}; assert(isequal(size(Xi), [C T]), 'imgs4d: size mismatch');
    Im4d(:,:,1,i) = single(Xi);
end
Y1 = Y;
end

function [W, mC, expl, K] = pca_fit_channels_cell(Xcells, varKeep)
% PCA across channels: concat all time samples from TRAIN
C = size(Xcells{1},1);
Xcat = []; Xcat = cat(2, Xcells{:});   % C × sum(T)
mC = mean(Xcat, 2);
X0 = Xcat - mC;
[U,S,~] = svd(double(X0), 'econ');
s = diag(S);
expl = 100*(s.^2)/sum(s.^2);
K = find(cumsum(expl)>=varKeep, 1, 'first'); if isempty(K), K=size(U,2); end
W = U;
end

function Zcells = pca_apply_channels_cell(Xcells, W, mC, K)
Zcells = cell(size(Xcells));
Wk = W(:,1:K);
for i=1:numel(Xcells)
    x = single(Xcells{i});
    Zcells{i} = Wk'*(x - mC);  % K × T
end
end

function [alpha, BA, cm] = pick_fusion_alpha(net1, net2, Xva1, ImVa, YVal, cats)
[~, s1] = classify(net1, Xva1); p1 = s1(:, strcmp(cats,'Spike'));
[~, s2] = classify(net2, ImVa); p2 = s2(:, strcmp(cats,'Spike'));
alphas = linspace(0,1,41); bestBA=-Inf; bestA=0.5; bestCM=[];
for a = alphas
    p = a*p1 + (1-a)*p2;
    Yp = categorical( (p>=0.5)+1, [1 2], cats);
    [acc, BAcur, cmcur] = metrics(YVal,Yp);
    if BAcur > bestBA, bestBA=BAcur; bestA=a; bestCM=cmcur; end
end
alpha = bestA; BA = bestBA; cm = bestCM;
end

function [acc, BA, cm] = metrics(Ytrue, Yhat)
cats = categories(Ytrue); Ytrue = reordercats(Ytrue,cats);
if isa(Yhat,'categorical'), Yhat = reordercats(Yhat,cats);
else, Yhat = categorical(Yhat, cats, 'Ordinal', false); end
cm  = confusionmat(Ytrue, Yhat, 'Order', cats);
acc = sum(diag(cm))/sum(cm,'all');
rec = diag(cm)./max(1,sum(cm,2));
BA  = mean(rec);
end

function [mu,sd] = fit_channel_stats(C)
ch = size(C{1},1); sumC = zeros(ch,1); sumSqC = zeros(ch,1); nTot=0;
for i=1:numel(C)
    xi=double(C{i}); sumC=sumC+sum(xi,2); sumSqC=sumSqC+sum(xi.^2,2); nTot=nTot+size(xi,2);
end
mu = sumC./nTot; varC = (sumSqC./nTot) - mu.^2; sd = sqrt(max(varC,1e-12));
end

function [med,madv] = fit_channel_stats_robust(C)
ch = size(C{1},1); med = zeros(ch,1); madv = ones(ch,1);
for r = 1:ch
    rows = cellfun(@(z) double(z(r,:)), C, 'UniformOutput', false);
    rows = rows(~cellfun(@isempty, rows));
    if isempty(rows), med(r)=0; madv(r)=1; continue; end
    x = [rows{:}];
    m = median(x,'omitnan');
    med(r)=m; madv(r)=max(1e-12, 1.4826*median(abs(x-m),'omitnan'));
end
end

function [Xc,Yc,Gc] = cap_by_group(X,Y,G,capPerClass)
uG = unique(G,'stable'); Xc={}; Yc=categorical; Gc={};
for gi = 1:numel(uG)
    idxG = find(G==uG{gi});
    for c = categories(Y).'
        idx = idxG(Y(idxG)==c);
        if isempty(idx), continue; end
        if numel(idx) > capPerClass, idx = idx(randperm(numel(idx), capPerClass)); end
        Xc = [Xc; X(idx)]; Yc = [Yc; Y(idx)]; Gc = [Gc; G(idx)]; %#ok<AGROW>
    end
end
Yc = categorical(Yc, categories(Y), 'Ordinal', false);
end

function lgraph = resnet1d_layers(inputSize, numClasses)
lgraph = layerGraph();
lgraph = addLayers(lgraph, sequenceInputLayer(inputSize,'Normalization','none','Name','in'));
stem = [ convolution1dLayer(7,32,'Padding','same','Stride',1,'Name','stem_conv','WeightsInitializer','he')
    batchNormalizationLayer('Name','stem_bn')
    reluLayer('Name','stem_relu') ];
lgraph = addLayers(lgraph, stem);
lgraph = connectLayers(lgraph,'in','stem_conv');
prevName='stem_relu'; prevCh=32;
cfg = [64 1;64 2;64 4; 128 1;128 2;128 4; 256 1;256 2;256 4];
for i=1:size(cfg,1)
    [lgraph, prevName, prevCh] = addResBlock1d(lgraph, prevName, prevCh, cfg(i,1), 3, cfg(i,2), sprintf('b%d',i), 0.1);
end
head = [ globalAveragePooling1dLayer("Name","gap")
    dropoutLayer(0.2,"Name","head_drop")
    fullyConnectedLayer(numClasses,"Name","fc")
    softmaxLayer("Name","sm")
    classificationLayer("Name","cls") ];
lgraph = addLayers(lgraph, head);
lgraph = connectLayers(lgraph, prevName,'gap');
end

function [lgraph, outName, outCh] = addResBlock1d(lgraph, inName, inCh, filters, k, dil, blk, dropP)
main = [ convolution1dLayer(k, filters, 'Padding','same','DilationFactor',dil,'Name',[blk '_conv1'],'WeightsInitializer','he')
    batchNormalizationLayer('Name',[blk '_bn1'])
    reluLayer('Name',[blk '_relu1'])
    dropoutLayer(dropP,'Name',[blk '_drop'])
    convolution1dLayer(k, filters, 'Padding','same','DilationFactor',1,'Name',[blk '_conv2'],'WeightsInitializer','he')
    batchNormalizationLayer('Name',[blk '_bn2']) ];
lgraph = addLayers(lgraph, main);
if inCh ~= filters
    skip = [ convolution1dLayer(1, filters, 'Padding','same','Name',[blk '_proj'],'WeightsInitializer','he')
        batchNormalizationLayer('Name',[blk '_bnskip']) ];
    lgraph = addLayers(lgraph, skip);
    skipOut = [blk '_bnskip']; lgraph = connectLayers(lgraph, inName, [blk '_proj']);
else
    skipOut = inName;
end
addL  = additionLayer(2,'Name',[blk '_add']);
relu2 = reluLayer('Name',[blk '_relu2']);
lgraph = addLayers(lgraph, addL); lgraph = addLayers(lgraph, relu2);
lgraph = connectLayers(lgraph, inName,       [blk '_conv1']);
lgraph = connectLayers(lgraph, [blk '_bn2'], [blk '_add/in1']);
lgraph = connectLayers(lgraph, skipOut,      [blk '_add/in2']);
lgraph = connectLayers(lgraph, [blk '_add'], [blk '_relu2']);
outName = [blk '_relu2']; outCh = filters;
end

function [Im4d, Y1] = imgs4d_depth(Xcells, Y)
% Xcells: {N} each C×T×D  -> Im4d: C×T×D×N
N = numel(Xcells);
[C,T,D] = size(Xcells{1});
Im4d = zeros(C,T,D,N,'single');
for i = 1:N
    Xi = Xcells{i};
    assert(ndims(Xi)==3 && isequal(size(Xi),[C,T,D]), 'imgs4d_depth: size mismatch');
    Im4d(:,:,:,i) = single(Xi);
end
Y1 = Y;
end

