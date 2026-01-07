%% ------------------------------------------------------------------------
% Build Spike/NoSpike dataset + leak-free split + normalization + baseline CNN
% -------------------------------------------------------------------------
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/function')

root   = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files';
raters = {'Adi','Josh','Manoj','Pradeep'};

% -------------------- User settings --------------------------------------
cropSamp  = 1:306;     % keep first 306 sensors (your current choice)
useSingle = true;      % store X as single to save RAM
Fs_default = 128;      % fallback if Fs not stored
rng(1);                % reproducible split

% -------------------- Collect all trials ---------------------------------
allX = [];
ally = [];
subj_id = {};
rater_id = {};
base_id  = {};
Fs_list  = [];
pairs = {};

for r = 1:numel(raters)
    d = fullfile(root, raters{r});
    spk = dir(fullfile(d, '*_spike_events.mat'));

    for k = 1:numel(spk)
        spk_fn = fullfile(spk(k).folder, spk(k).name);

        base = erase(spk(k).name, '_spike_events.mat');
        nos_fn = fullfile(d, [base '_nospike_events.mat']);

        if ~isfile(nos_fn)
            fprintf('Missing nospike pair for: %s\n', spk(k).name);
            continue;
        end

        [Es, metas] = load_epochs_mat(spk_fn);
        [En, metan] = load_epochs_mat(nos_fn);

        % --- Expect: trials × time × chan ---
        % Crop time safely (in case a file is shorter)
        Es = Es(:, cropSamp(cropSamp <= size(Es,2)), :);
        En = En(:, cropSamp(cropSamp <= size(En,2)), :);

        % Convert to: trials × chan × time  (N×C×T)
        Xs = permute(Es, [1 3 2]);
        Xn = permute(En, [1 3 2]);

        % Labels
        y = [ones(size(Xs,1),1); zeros(size(Xn,1),1)];

        % Concatenate
        X = cat(1, Xs, Xn);
        if useSingle, X = single(X); end

        allX = cat(1, allX, X);
        ally = cat(1, ally, y);

        % IDs
        sid = extract_subject_id(base);  % subject-level grouping (prevents leakage)
        subj_id = [subj_id; repmat({sid}, numel(y), 1)];
        rater_id = [rater_id; repmat(raters(r), numel(y), 1)];
        base_id  = [base_id;  repmat({base}, numel(y), 1)];

        % Fs bookkeeping
        Fs_s = getFs(metas, Fs_default);
        Fs_n = getFs(metan, Fs_default);
        Fs_list = [Fs_list; repmat(Fs_s, numel(y), 1)]; %#ok<AGROW>
        if Fs_s ~= Fs_n
            warning('Fs mismatch for %s (%g vs %g)', base, Fs_s, Fs_n);
        end

        pairs(end+1,:) = {raters{r}, base, spk_fn, nos_fn, ...
            size(Xs,1), size(Xn,1), size(Xs,2), size(Xs,3)}; 
    end
end

T = cell2table(pairs, 'VariableNames', ...
    {'Rater','Base','SpikeFile','NoSpikeFile','nSpike','nNoSpike','nChan','nTime'});
disp(T)

fprintf('Final dataset: X = [%d trials × %d chan × %d time], y = [%d×1]\n', ...
    size(allX,1), size(allX,2), size(allX,3), numel(ally));

%% -------------------- Quick sanity plots --------------------------------
Fs = mode(Fs_list(~isnan(Fs_list)));
if isempty(Fs), Fs = Fs_default; end
t = (0:size(allX,3)-1)/Fs;

trial = 1;
figure; plot(t, squeeze(allX(trial,1,:)));
xlabel('Time (s)'); ylabel('Amplitude');
title(sprintf('Trial %d, Channel 1, y=%d', trial, ally(trial)));
grid on;

figure; imagesc(t, 1:size(allX,2), squeeze(allX(trial,:,:)));
axis tight; xlabel('Time (s)'); ylabel('Channel');
title('Trial heatmap (channels × time)'); colorbar;

% Mean template (spike vs nospike), helps confirm signal differences
Xs_all = allX(ally==1,:,:);
Xn_all = allX(ally==0,:,:);
mS = squeeze(mean(Xs_all,1,'omitnan'));  % C×T
mN = squeeze(mean(Xn_all,1,'omitnan'));  % C×T
figure; imagesc(t, 1:size(allX,2), mS); axis tight; colorbar;
title('Mean Spike (channels × time)'); xlabel('Time (s)'); ylabel('Channel');
figure; imagesc(t, 1:size(allX,2), mN); axis tight; colorbar;
title('Mean NoSpike (channels × time)'); xlabel('Time (s)'); ylabel('Channel');

%% -------------------- Leak-free split by SUBJECT ------------------------
subjects = unique(subj_id);
nSubj = numel(subjects);

% 70/15/15 split
idx = randperm(nSubj);
nTrain = round(0.70*nSubj);
nVal   = round(0.15*nSubj);
trainSubj = subjects(idx(1:nTrain));
valSubj   = subjects(idx(nTrain+1:nTrain+nVal));
testSubj  = subjects(idx(nTrain+nVal+1:end));

isTrain = ismember(subj_id, trainSubj);
isVal   = ismember(subj_id, valSubj);
isTest  = ismember(subj_id, testSubj);

fprintf('Subjects: train=%d, val=%d, test=%d\n', numel(trainSubj), numel(valSubj), numel(testSubj));
fprintf('Trials:    train=%d, val=%d, test=%d\n', sum(isTrain), sum(isVal), sum(isTest));

%% -------------------- Normalize using TRAIN stats only -------------------
[X_train_norm, mu, sig] = normalize_by_train(allX(isTrain,:,:));
X_val_norm  = apply_norm(allX(isVal,:,:),  mu, sig);
X_test_norm = apply_norm(allX(isTest,:,:), mu, sig);

y_train = ally(isTrain);
y_val   = ally(isVal);
y_test  = ally(isTest);

%% -------------------- Save DL-ready dataset ------------------------------
dataset = struct();
dataset.X_train = X_train_norm; dataset.y_train = y_train;
dataset.X_val   = X_val_norm;   dataset.y_val   = y_val;
dataset.X_test  = X_test_norm;  dataset.y_test  = y_test;

dataset.subject_id_train = subj_id(isTrain);
dataset.subject_id_val   = subj_id(isVal);
dataset.subject_id_test  = subj_id(isTest);

dataset.rater_train = rater_id(isTrain);
dataset.rater_val   = rater_id(isVal);
dataset.rater_test  = rater_id(isTest);

dataset.base_train = base_id(isTrain);
dataset.base_val   = base_id(isVal);
dataset.base_test  = base_id(isTest);

dataset.mu = mu; dataset.sig = sig;
dataset.Fs = Fs;

save('dataset_DL_ready.mat', 'dataset', '-v7.3');
disp('Saved: dataset_DL_ready.mat');

%% -------------------- OPTIONAL: baseline CNN in MATLAB -------------------
% Treat each trial as an "image": [channels × time × 1]
% X is currently N×C×T -> need C×T×1×N
Xtr4 = to4D(dataset.X_train);  Ytr = toCat(dataset.y_train);
Xva4 = to4D(dataset.X_val);    Yva = toCat(dataset.y_val);

C = size(dataset.X_train,2);
Tt = size(dataset.X_train,3);

% Class weights (handle imbalance)
nSpike = sum(dataset.y_train==1);
nNo   = sum(dataset.y_train==0);
wSpike = (nSpike+nNo)/(2*nSpike+eps);
wNo    = (nSpike+nNo)/(2*nNo+eps);
classWeights = [wNo, wSpike];  % order matches categorical below: NoSpike, Spike

layers = [
    imageInputLayer([C Tt 1], 'Normalization','none','Name','in')

    convolution2dLayer([5 9], 16, 'Padding','same','Name','conv1')
    batchNormalizationLayer('Name','bn1')
    reluLayer('Name','relu1')
    maxPooling2dLayer([2 2], 'Stride',[2 2], 'Name','pool1')

    convolution2dLayer([3 7], 32, 'Padding','same','Name','conv2')
    batchNormalizationLayer('Name','bn2')
    reluLayer('Name','relu2')
    maxPooling2dLayer([2 2], 'Stride',[2 2], 'Name','pool2')

    convolution2dLayer([3 5], 64, 'Padding','same','Name','conv3')
    batchNormalizationLayer('Name','bn3')
    reluLayer('Name','relu3')

    globalAveragePooling2dLayer('Name','gap')
    dropoutLayer(0.2,'Name','drop')
    fullyConnectedLayer(2,'Name','fc')

    softmaxLayer('Name','sm')
    classificationLayer('Name','cls', 'Classes', categorical({'NoSpike','Spike'}), 'ClassWeights', classWeights)
];

opts = trainingOptions('adam', ...
    'MaxEpochs', 25, ...
    'MiniBatchSize', 64, ...
    'Shuffle','every-epoch', ...
    'ValidationData', {Xva4, Yva}, ...
    'ValidationFrequency', 50, ...
    'Plots','training-progress', ...
    'Verbose', true);

net = trainNetwork(Xtr4, Ytr, layers, opts);

% Evaluate on test
Xte4 = to4D(dataset.X_test);
Yte  = toCat(dataset.y_test);
Yhat = classify(net, Xte4);
acc = mean(Yhat == Yte);
fprintf('Test accuracy: %.3f\n', acc);

%% -------------------- Local helper functions ----------------------------
function Fs = getFs(meta, Fs_default)
    Fs = Fs_default;
    if isstruct(meta) && isfield(meta,'Fs') && ~isempty(meta.Fs)
        Fs = double(meta.Fs);
        if numel(Fs) > 1, Fs = Fs(1); end
    end
end

function sid = extract_subject_id(base)
    % base like: alioto_courtney_Run04_spont_supine_raw_t_sss_ecgClean
    tok = regexp(base, '^(.*)_Run\d+', 'tokens', 'once');
    if ~isempty(tok)
        sid = tok{1};
    else
        sid = base; % fallback: group by full base if pattern doesn't match
    end
end

function [Xn, mu, sig] = normalize_by_train(X)
    % X: N×C×T
    [N,C,T] = size(X); %#ok<ASGLU>
    tmp = permute(X, [2 1 3]);     % C×N×T
    tmp = reshape(tmp, C, []);     % C×(N*T)
    mu  = mean(tmp, 2, 'omitnan');
    sig = std(tmp, 0, 2, 'omitnan');
    sig(sig==0) = 1;

    Xn = apply_norm(X, mu, sig);
end

function Xn = apply_norm(X, mu, sig)
    % X: N×C×T, mu/sig: C×1
    mu3  = reshape(mu,  1, [], 1);
    sig3 = reshape(sig, 1, [], 1);
    Xn = (X - mu3) ./ sig3;
end

function X4 = to4D(X)
    % X: N×C×T -> C×T×1×N
    X4 = permute(X, [2 3 1]);          % C×T×N
    X4 = reshape(X4, size(X4,1), size(X4,2), 1, size(X4,3));
end

function Y = toCat(y)
    % y: 0/1 -> categorical
    Y = categorical(y, [0 1], {'NoSpike','Spike'});
end