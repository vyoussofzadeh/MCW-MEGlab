function run_spike_classifier2(varargin)
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

% Optional skip list
isSkip = @(f) any(cellfun(@(s) contains(f,s), opts.skipList));
filesS = filesS(~arrayfun(@(f) isSkip(f.name), filesS));
filesN = filesN(~arrayfun(@(f) isSkip(f.name), filesN));

%% ----------------------- BUILD TRIAL CELL ARRAYS ------------------------
% Returns epoch cells and group IDs (subject/session) for leakage-safe split
proc_class = @(flg,files) build_cells(files,flg,@extract_group_id);
[dataSpike, groupsS]  = proc_class(1, filesS);
[dataNoSpike, groupsN]= proc_class(0, filesN);

% remove NaNs class-wise (keep group alignment)
clean  = @(C,G) deal(C(~cellfun(@(x) any(isnan(x(:))),C)), G(~cellfun(@(x) any(isnan(x(:))),C)));
[Dspk, groupsS]   = clean(dataSpike, groupsS);
[Dnspk, groupsN]  = clean(dataNoSpike, groupsN);

data   = [Dnspk; Dspk];
labels = categorical([zeros(numel(Dnspk),1); ones(numel(Dspk),1)], [0 1], {'NoSpike','Spike'});
groups = [groupsN; groupsS];

%%
% save('DD.mat','data','labels','groups')

%% ------------------------ GROUPED TRAIN / VAL SPLIT ---------------------
[trainMask, valMask] = group_holdout_split(groups, opts.valFrac, opts.seed);
XTrain = data(trainMask);   YTrain = labels(trainMask);
XVal   = data(valMask);     YVal   = labels(valMask);

%%
for i=1:length(XTrain)
    XTrain{i} = XTrain{i}(:,1:67);
end

%% ----------------------- CHANNEL-WISE NORMALIZATION ---------------------
% Fit on training only; supports 'zscore' or 'robust' (median/MAD)
switch lower(opts.norm)
    case 'zscore'
        [mu, sd] = fit_channel_stats(XTrain);
        zscoref  = @(x) (x-mu)./sd;      % protect sd inside fitter
        normStats = struct('type','zscore','mu',mu,'sd',sd);
    case 'robust'
        [med, madv] = do_fit_channel_stats_robust(XTrain);
        zscoref     = @(x) (x-med)./max(madv,1e-12);
        normStats = struct('type','robust','med',med,'madv',madv);
        mu = []; sd = []; % for backward compatibility variables
    otherwise
        error('Unknown opts.norm: %s', opts.norm);
end
XTrain2 = cellfun(zscoref, XTrain,'UniformOutput',false);
XVal2   = cellfun(zscoref, XVal  ,'UniformOutput',false);

% Optional light augmentation (training only)
if opts.augment
    XTrain2 = cellfun(@augment_epoch, XTrain2,'UniformOutput',false);
end

% Append positional channels for Transformer only (no toolbox deps)
if strcmpi(opts.arch,'transformer')
    addPos = @(x) [x; 0.1*make_posenc(size(x,2),32)]; % scaled PE (32 extra ch)
    XTrain2 = cellfun(addPos,XTrain2,'UniformOutput',false);
    XVal2   = cellfun(addPos,XVal2  ,'UniformOutput',false);
end

% Balance classes by oversampling minority in training only
[XTrain2, YTrain] = balance_oversample(XTrain2, YTrain);

%% ------------------------- BUILD MODEL & OPTIONS ------------------------
% Sanity: consistent channel counts
chCounts = cellfun(@(x) size(x,1), XTrain2);
assert(~isempty(chCounts),'No training epochs left after cleaning.');
uCh = unique(chCounts);
assert(numel(uCh)==1, 'Inconsistent channel counts across epochs: %s', mat2str(uCh));
inputSize = uCh(1);

padArgs = {"SequenceLength","longest"};   % keep all samples; pad per mini-batch

%%

% apply to train/val AFTER your z-scoring & (optional) augmentation:
% XTrain2 = cellfun(@add_global_feats, XTrain2, 'UniformOutput', false);
% XVal2   = cellfun(@add_global_feats, XVal2,   'UniformOutput', false);
% 
% % and recompute inputSize accordingly before defining layers:
% inputSize = size(XTrain2{1},1);

%%
switch lower(opts.arch)
    case 'lstm'
mb = 128;  % try 64 if you still overfit
valFreq = max(5, ceil(numel(XTrain2)/mb));

trainOpts = trainingOptions('adam', ...
    'InitialLearnRate',3e-4, ...
    'MaxEpochs',120, ...
    'MiniBatchSize',mb, ...
    'L2Regularization',1e-3, ...
    'GradientThresholdMethod','l2norm', ...
    'GradientThreshold',1, ...
    'SequenceLength','longest', ...
    'SequencePaddingValue',0, ...
    'SequencePaddingDirection','right', ...
    'Shuffle','every-epoch', ...
    'ValidationData',{XVal2,YVal}, ...
    'ValidationFrequency',valFreq, ...
    'ValidationPatience',12, ...
    'ExecutionEnvironment', ternary(opts.useGPU,'gpu','cpu'), ...
    'Verbose',opts.verbose, ...
    'Plots','training-progress', ...
    'CheckpointPath',opts.modelDir);

layers = [
    sequenceInputLayer(inputSize,'Normalization','none')
    bilstmLayer(128,'OutputMode','sequence')
    dropoutLayer(0.3)                              % a touch more dropout
    bilstmLayer(64,'OutputMode','sequence')        % <-- was 'last'
    layerNormalizationLayer                        % stabilize activations
    globalMaxPooling1dLayer                        % <-- robust to jitter
    dropoutLayer(0.3)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer
];





        % mb = 64; lr = 1e-3; maxE = 120;
        % valFreq = max(5, ceil(numel(XTrain2)/mb));
        % layers = [
        %     sequenceInputLayer(inputSize,'Normalization','none')
        %     bilstmLayer(128,'OutputMode','sequence')
        %     dropoutLayer(0.2)
        %     bilstmLayer(64,'OutputMode','last')
        %     dropoutLayer(0.2)
        %     fullyConnectedLayer(2)
        %     softmaxLayer
        %     classificationLayer];
        % trainOpts = trainingOptions('adam', ...
        %     'InitialLearnRate',lr, ...
        %     'MaxEpochs',maxE, ...
        %     'MiniBatchSize',mb, ...
        %     'GradientThresholdMethod','l2norm', ...
        %     'GradientThreshold',1, ...
        %     'Shuffle','every-epoch', ...
        %     padArgs{:}, ...
        %     'ValidationData',{XVal2,YVal}, ...
        %     'ValidationFrequency',valFreq, ...
        %     'ValidationPatience',15, ...
        %     'ExecutionEnvironment', ternary(opts.useGPU,'gpu','cpu'), ...
        %     'Verbose',opts.verbose, ...
        %     'Plots','training-progress', ...
        %     'CheckpointPath',opts.modelDir);
        net = trainNetwork(XTrain2,YTrain,layers,trainOpts);
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
        net = trainNetwork(XTrain2,YTrain,lg,trainOpts);
        mdlName = 'grutd_spike_classifier.mat';

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
        net = trainNetwork(XTrain2,YTrain,lg,trainOpts);
        mdlName = 'transformer_spike_classifier.mat';

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
meta = struct('matlabVersion',version,'gitHash',git_hash,'timestamp',datestr(now,'yyyy-mm-dd HH:MM:SS'));

save(fullfile(opts.modelDir,mdlName), 'net','normStats','mu','sd','opts','trainOpts','cm','acc','balAcc','f1','auc','meta');
save(fullfile(opts.modelDir,'run_meta.mat'),'meta');

% Figures
if opts.savePlots
    fig1=figure; plotconfusion(YVal,YP); title(sprintf('Val Acc=%.2f, BA=%.2f, AUC=%.2f',acc,balAcc,auc));
    saveas(fig1, fullfile(opts.modelDir,'confusion_val.png')); close(fig1);
    if ~isempty(rocX)
        fig2=figure; plot(rocX,rocY); xlabel('FPR'); ylabel('TPR'); title(sprintf('ROC (AUC=%.3f)',auc)); grid on;
        saveas(fig2, fullfile(opts.modelDir,'roc_val.png')); close(fig2);
    end
    if ~isempty(XTrain2)
        fig3=figure('Name','Example epoch'); clf
        plot(XTrain2{1}'); title('First normalized epoch'); xlabel('Time');
        saveas(fig3, fullfile(opts.modelDir,'example_epoch.png')); close(fig3);
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

function [trainMask,valMask] = group_holdout_split(groups, valFrac, seed)
% Splits by unique group IDs to prevent leakage between train/val
if nargin<3, seed = 0; end
rng(seed,'twister');
u = unique(groups,'stable');
nVal = max(1, round(valFrac * numel(u)));
idx = randperm(numel(u));
valGroups = u(idx(1:nVal));
valMask = ismember(groups, valGroups);
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
% Robust per-channel stats (median, MAD->sigma)
ch = size(C{1},1); med = zeros(ch,1); madv = zeros(ch,1);
for r=1:ch
    x = cell2mat(cellfun(@(z) double(z(r,:)), C,'UniformOutput',false));
    med(r)  = median(x,'omitnan');
    madv(r) = 1.4826*mad(x,1);  % MAD to sigma
end
madv = max(madv,1e-12);
end

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

function Xaug = augment_epoch(x)
% Very light augmentation to improve generalization
if rand<0.5, x = x + 0.005*randn(size(x)); end        % Gaussian jitter
if rand<0.3                                          % time mask
    T = size(x,2); w = randi([max(1,round(0.02*T)), max(1,round(0.08*T))]);
    s = randi([1, max(1, T-w+1)]); x(:, s:min(T,s+w-1)) = 0;
end
Xaug = x;
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
% Extended parser with arch/norm/augment/seed
p = inputParser; p.FunctionName = 'run_spike_classifier';
p.CaseSensitive = false; p.PartialMatching = true; p.KeepUnmatched = false;
addParameter(p,'ftRoot',def.ftRoot);
addParameter(p,'codeRoot',def.codeRoot);
addParameter(p,'spikeDir',def.spikeDir);
addParameter(p,'noSpikeDir',def.noSpikeDir);
addParameter(p,'skipList',def.skipList);
addParameter(p,'useGPU',def.useGPU);
addParameter(p,'modelDir',def.modelDir);
addParameter(p,'savePlots',def.savePlots);
addParameter(p,'verbose',def.verbose);
addParameter(p,'valFrac',def.valFrac);
validArch = @(s) any(strcmpi(s,{'lstm','gru_td','transformer'}));
addParameter(p,'arch',def.arch, validArch);
addParameter(p,'dModel',def.dModel, @(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
addParameter(p,'numHeads',def.numHeads, @(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
addParameter(p,'norm',def.norm, @(s) any(strcmpi(s,{'zscore','robust'})));
addParameter(p,'augment',def.augment, @(x) islogical(x) || isnumeric(x));
addParameter(p,'seed',def.seed, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:}); out = p.Results;
if strcmpi(out.arch,'transformer')
    assert(mod(out.dModel,out.numHeads)==0, 'dModel (%d) must be divisible by numHeads (%d).', out.dModel, out.numHeads);
end
rng(out.seed,'twister'); % reset seed if user overrode
end

function setup_environment(o)
addpath(o.ftRoot); ft_defaults;
addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/functions_new')));
addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/helper')));
end

function r = ternary(c,a,b); if c, r=a; else, r=b; end; end


function [med,madv] = do_fit_channel_stats_robust(C)
% Robust per-channel stats (median, MADsigma) for variable-length epochs
% Concatenate each channel as a single long row vector across epochs.
if isempty(C) || isempty(C{1})
    error('fit_channel_stats_robust:EmptyTraining','XTrain is empty; cannot fit robust stats.');
end
ch = size(C{1},1); 
med = zeros(ch,1); 
madv = zeros(ch,1);

for r = 1:ch
    rows = cellfun(@(z) double(z(r,:)), C, 'UniformOutput', false);
    % force row-vector shape and drop empties
    rows = cellfun(@(v) reshape(v,1,[]), rows, 'UniformOutput', false);
    rows = rows(~cellfun(@isempty, rows));

    if isempty(rows)
        med(r)  = 0; 
        madv(r) = 1;
        continue;
    end

    x = [rows{:}];  % 1 x (sum T)
    if isempty(x) || all(isnan(x))
        med(r)  = 0; 
        madv(r) = 1;
    else
        m = median(x,'omitnan');
        med(r)  = m;
        madv(r) = 1.4826 * median(abs(x - m),'omitnan'); % MAD->s
    end
end

madv = max(madv,1e-12);  % protect against divide-by-zero
end

% function Xaug = augment_epoch(x)
%     if rand<0.5, x = x + 0.005*randn(size(x)); end   % jitter
%     if rand<0.3
%         T = size(x,2); w = randi([round(0.02*T), round(0.08*T)]);
%         s = randi([1, T-w+1]); x(:, s:s+w-1) = 0;     % time mask
%     end
%     Xaug = x;
% end

function y = add_global_feats(x)
    gfp    = sqrt(mean(x.^2,1));
    absmax = max(abs(x),[],1);
    energy = sum(x.^2,1);
    y = [x; gfp; absmax; energy];
end

