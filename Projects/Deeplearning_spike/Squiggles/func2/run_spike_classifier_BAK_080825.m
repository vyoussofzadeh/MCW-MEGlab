%% Spike-Detection MEG Pipeline  Spike vs Non-Spike Classifier (July 22 2025)
% Deep-learning workflow to classify inter-ictal MEG epochs as **Spike** or **NoSpike**.
% Author: MCW MEG Lab  V. Youssofzadeh <vyoussofzadeh@mcw.edu>
% -------------------------------------------------------------------------
% Usage:
%   1. Edit the USER SETTINGS section.
%   2. Run the script or call `run_spike_classifier` from MATLAB.
% -------------------------------------------------------------------------
function run_spike_classifier(varargin)

%% --------------------------- USER SETTINGS -----------------------------
opts.ftRoot      = '/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
opts.codeRoot    = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git';
opts.spikeDir    = '/data/MEG/Research/SpikeDectection/Epil_annotated_data/annotated_data';
opts.noSpikeDir  = '/data/MEG/Research/SpikeDectection/Epil_annotated_data/annotated_data_nospike';
opts.skipList    = {};                 % e.g. {'pilot','bad_subj'}
opts.useGPU      = false;               % auto-downgrades if incompatible
opts.modelDir    = fullfile(opts.spikeDir,'models_classifier');
opts.savePlots   = true;
opts.verbose     = true;
opts.valFrac     = 0.15;
opts.arch        = 'transformer';     % <--- 'lstm' | 'gru_td' | 'transformer'
opts.dModel      = 128;        % transformer width (ignored by LSTM/GRU)
opts.numHeads    = 8;          % transformer heads (ignored by LSTM/GRU)
rng(0,'twister');

% Allow name-value overrides
opts = parse_inputs(opts,varargin{:});

%% ----------------------------- ENV SETUP -------------------------------
setup_environment(opts);

gpuOK = false;
if opts.useGPU
    try g=gpuDevice; gpuOK = g.ComputeCapability >= "3.0"; catch, end
end
if ~gpuOK, opts.useGPU = false; disp('Training on CPU (GPU unavailable/old).'); end
if ~exist(opts.modelDir,'dir'); mkdir(opts.modelDir); end

%% ------------------------- LOAD MAT FILE LISTS -------------------------
add_files  = @(p) rdir(fullfile(p,'*.mat'));
filesS = add_files(opts.spikeDir);
filesN = add_files(opts.noSpikeDir);
assert(~isempty(filesS) && ~isempty(filesN),'No .mat files found.');

%% Optional skip list
isSkip = @(f) any(cellfun(@(s) contains(f,s), opts.skipList));
filesS = filesS(~arrayfun(@(f) isSkip(f.name), filesS));
filesN = filesN(~arrayfun(@(f) isSkip(f.name), filesN));

%% ----------------------- BUILD TRIAL CELL ARRAYS -----------------------
refLabels   = [];
proc_class  = @(flg,files) build_cells(files,flg);
[dataSpike, ~]   = proc_class(1, filesS);
[dataNoSpike, ~]         = proc_class(0, filesN);

% remove NaNs class-wise
clean  = @(C) C(~cellfun(@(x) any(isnan(x(:))),C));
Dspk   = clean(dataSpike);  Dnspk = clean(dataNoSpike);

labels = categorical([zeros(numel(Dnspk),1); ones(numel(Dspk),1)], ...
    [0 1], {'NoSpike','Spike'});
data   = [Dnspk; Dspk];

%% ------------------------ TRAIN / VAL SPLIT ----------------------------
cv = cvpartition(numel(labels),'HoldOut',opts.valFrac);

XTrain = data(training(cv));   YTrain = labels(training(cv));
XVal   = data(test(cv));       YVal   = labels(test(cv));

% -------- normalise (fit mu,sd on training only) --------
% Some epochs have 67-102 samples ? direct cat may fail. Compute
% channel-wise µ,s incrementally to avoid dimension-mismatch.
ch = size(XTrain{1},1);
sumC = zeros(ch,1);  sumSqC = zeros(ch,1);  nTot = 0;
for i = 1:numel(XTrain)
    xi      = XTrain{i};
    sumC    = sumC  + sum(xi,2);
    sumSqC  = sumSqC+ sum(xi.^2,2);
    nTot    = nTot  + size(xi,2);
end
mu = sumC ./ nTot;
varC = (sumSqC ./ nTot) - mu.^2;
sd = sqrt(max(varC, 1e-12));   % protect zero-variance chans
zscoref  = @(x) (x-mu)./sd;

XTrain2 = []; for i=1:length(XTrain), XTrain2{i} = XTrain{i}; end
XVal2 = []; for i=1:length(XVal), XVal2{i} = XVal{i}; end

XTrain2 = cellfun(zscoref, XTrain2,'UniformOutput',false);
XVal2   = cellfun(zscoref, XVal2  ,'UniformOutput',false);

% ---- Append positional channels for Transformer only (no toolbox deps) ----
if strcmpi(opts.arch,'transformer')
    addPos = @(x) [x; make_posenc(size(x,2),32)]; % 32 extra channels (sin/cos)
    XTrain2 = cellfun(addPos,XTrain2,'UniformOutput',false);
    XVal2   = cellfun(addPos,XVal2,'UniformOutput',false);
end

%% ------------------------- BUILD MODEL & OPTIONS -----------------------
inputSize = size(XTrain2{1},1);

switch lower(opts.arch)
    case 'lstm'
        layers = [
            sequenceInputLayer(inputSize,'Normalization','none')
            bilstmLayer(128,'OutputMode','sequence')
            dropoutLayer(0.2)
            bilstmLayer(64,'OutputMode','last')
            dropoutLayer(0.2)
            fullyConnectedLayer(2)
            softmaxLayer
            classificationLayer];
        
        trainOpts = trainingOptions('adam', ...
            'MaxEpochs',120, ...
            'MiniBatchSize',256, ...
            'InitialLearnRate',1e-3, ...
            'LearnRateSchedule','piecewise', ...
            'LearnRateDropFactor',0.5, ...
            'LearnRateDropPeriod',50, ...
            'GradientThreshold',1, ...
            'Shuffle','every-epoch', ...
            'ValidationData',{XVal2,YVal}, ...
            'ValidationFrequency',max(5,floor(numel(XTrain2)/256)*3), ...
            'ValidationPatience',10, ...
            'ExecutionEnvironment', ternary(opts.useGPU,'gpu','cpu'), ...
            'Verbose',opts.verbose, ...
            'Plots','training-progress', ...
            'CheckpointPath',opts.modelDir);
        
        net = trainNetwork(XTrain2,YTrain,layers,trainOpts);
        mdlName = 'lstm_spike_classifier.mat';
        
    case 'gru_td'  % GRU with time-distributed dense via fold/unfold (R2023b-safe)
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
        lg = connectLayers(lg,"fold/miniBatchSize","unfold/miniBatchSize"); % critical bookkeeping
        
        trainOpts = trainingOptions("adam", ...
            "MaxEpochs",100, ...
            "MiniBatchSize",64, ...
            "InitialLearnRate",1e-3, ...
            "ValidationData",{XVal2,YVal}, ...
            "ValidationFrequency",25, ...
            "ValidationPatience",10, ...
            "Shuffle","every-epoch", ...
            "ExecutionEnvironment", ternary(opts.useGPU,'gpu','cpu'), ...
            "Verbose",opts.verbose, ...
            "Plots","training-progress", ...
            "CheckpointPath",opts.modelDir);
        
        net = trainNetwork(XTrain2,YTrain,lg,trainOpts);
        mdlName = 'grutd_spike_classifier.mat';
        
    case 'transformer'
        
        
        % ---- (keep your pos-encoding block above this) ----
        % Sanity check: all epochs must have the same #channels
        chCounts = cellfun(@(x) size(x,1), XTrain2);
        assert(~isempty(chCounts),'No training epochs left after cleaning.');
        uCh = unique(chCounts);
        assert(numel(uCh)==1, ...
            'Inconsistent channel counts across epochs: %s', mat2str(uCh));
        
        inputSize = uCh(1);   % <--- now defined for all model branches
        
        
        
        dModel  = opts.dModel;
        numHeads= opts.numHeads;
        keyDim  = dModel/numHeads;  % must be integer
        
        lg = layerGraph();
        lg = addLayers(lg, sequenceInputLayer(inputSize,"Normalization","none","Name","in"));
        lg = addLayers(lg, fullyConnectedLayer(dModel,"Name","proj"));
        
        % --- Encoder block (pre-norm) ---
        lg = addLayers(lg, layerNormalizationLayer("Name","ln1"));
        lg = addLayers(lg, selfAttentionLayer(numHeads,keyDim, ...
            "DropoutProbability",0.1,"Name","mha1"));
        lg = addLayers(lg, dropoutLayer(0.1,"Name","drop_attn"));
        lg = addLayers(lg, additionLayer(2,"Name","add_attn"));
        
        lg = addLayers(lg, layerNormalizationLayer("Name","ln2"));
        lg = addLayers(lg, [
            fullyConnectedLayer(4*dModel,"Name","ff1")
            reluLayer("Name","relu1")
            dropoutLayer(0.1,"Name","drop_ff")
            fullyConnectedLayer(dModel,"Name","ff2")]);
        lg = addLayers(lg, additionLayer(2,"Name","add_ff"));
        
        % sequence ? label head
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
            "InitialLearnRate",3e-4, ...
            "MaxEpochs",150, ...
            "MiniBatchSize",128, ...
            "SequenceLength","shortest", ...  % trims per mini-batch
            "Shuffle","every-epoch", ...
            "ValidationData",{XVal2,YVal}, ...
            "ValidationFrequency",max(5,floor(numel(XTrain2)/128)), ...
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

%% ------------------------------ SAVE & EVAL ----------------------------
save(fullfile(opts.modelDir,mdlName), 'net','mu','sd','opts','trainOpts');

[YP,~] = classify(net,XVal2);
plotconfusion(YVal, YP);

%% ---------------------------- PLOTTING ---------------------------------
if opts.savePlots
    fig = figure('Name','Example epoch'); clf
    plot(zscoref(XTrain{1})'); title('First normalised epoch'); xlabel('Time');
    saveas(fig, fullfile(opts.modelDir,'example_epoch.png')); close(fig);
end
end % main

%% ------------------------- NETWORK ARCHITECTURE ------------------------
% inputSize = size(XTrain2{1},1);   % channels (rows)
% numHidden1 = 128;  numHidden2 = 64;
%
% layers = [
%     sequenceInputLayer(inputSize,'Normalization','none')
%     bilstmLayer(numHidden1,'OutputMode','sequence')
%     dropoutLayer(0.2)
%     bilstmLayer(numHidden2,'OutputMode','last')
%     dropoutLayer(0.2)
%     fullyConnectedLayer(2)
%     softmaxLayer
%     classificationLayer];
%
% %% model 1
% options = trainingOptions('adam', ...
%     'MaxEpochs',120, ...
%     'MiniBatchSize',256, ...
%     'InitialLearnRate',1e-3, ...
%     'LearnRateSchedule','piecewise', ...
%     'LearnRateDropFactor',0.5, ...
%     'LearnRateDropPeriod',50, ...
%     'GradientThreshold',1, ...
%     'Shuffle','every-epoch', ...
%     'ValidationData',{XVal2,YVal}, ...
%     'ValidationFrequency',floor(numel(XTrain2)/256)*3, ...
%     'ValidationPatience',10, ...
%     'ExecutionEnvironment',ternary(opts.useGPU,'gpu','cpu'), ...
%     'Verbose',opts.verbose, ...
%     'Plots','training-progress', ...
%     'CheckpointPath',opts.modelDir);
%
% %% model 2
% options = trainingOptions('adam', ...
%     'MaxEpochs',100, ...
%     'MiniBatchSize',64, ...
%     'InitialLearnRate',5e-3, ...
%     'LearnRateSchedule','piecewise', ...
%     'LearnRateDropPeriod',30, ...
%     'LearnRateDropFactor',0.5, ...
%     'Shuffle','every-epoch', ...
%     'ValidationData',{XVal2,YVal}, ...
%     'ValidationFrequency',25, ...
%     'ValidationPatience',10, ...
%     'ExecutionEnvironment', ternary(opts.useGPU,'gpu','cpu'), ...
%     'Verbose',opts.verbose, ...
%     'Plots','training-progress', ...
%     'CheckpointPath',opts.modelDir);
%
%
% %% ------------------------------ TRAIN ----------------------------------
% net = trainNetwork(XTrain2, YTrain, layers, options);
% save(fullfile(opts.modelDir,'lstm_spike_classifier.mat'), 'net','mu','sd','options');
%
% %%
% [YP,~] = classify(net,XVal2);
% plotconfusion(YVal, YP);

%% model 3
% %------------------------------------------------------
% % Build the layer graph
% %------------------------------------------------------
% inputSize = size(XTrain2{1},1);
%
% lg = layerGraph();
%
% % 1) input and fold
% lg = addLayers(lg, sequenceInputLayer(inputSize,...
%                       "Normalization","none","Name","input"));
% lg = addLayers(lg, sequenceFoldingLayer("Name","fold"));
%
% % 2) time-distributed dense (FC + ReLU)
% lg = addLayers(lg, fullyConnectedLayer(32,"Name","fc_td"));
% lg = addLayers(lg, reluLayer("Name","relu_td"));
%
% % 3) unfold to restore time axis
% lg = addLayers(lg, sequenceUnfoldingLayer("Name","unfold"));
% lg = addLayers(lg, flattenLayer("Name","flat"));   % tidy up dims
%
% % 4) recurrent pool + head
% lg = addLayers(lg, [
%         gruLayer(64,"OutputMode","last","Name","gru")
%         dropoutLayer(0.2,"Name","drop")
%         fullyConnectedLayer(2,"Name","fc")
%         softmaxLayer("Name","sm")
%         classificationLayer("Name","cls")]);
%
% % ---------- connect data path ----------
% lg = connectLayers(lg,"input","fold/in");
% lg = connectLayers(lg,"fold/out","fc_td");
% lg = connectLayers(lg,"fc_td","relu_td");
% lg = connectLayers(lg,"relu_td","unfold/in");
% lg = connectLayers(lg,"unfold/out","flat");
% lg = connectLayers(lg,"flat","gru");
%
% % ---------- critical bookkeeping path ----------
% lg = connectLayers(lg,"fold/miniBatchSize","unfold/miniBatchSize");
%
% %------------------------------------------------------
% % Train
% %------------------------------------------------------
% options = trainingOptions("adam", ...
%     "MaxEpochs",100, ...
%     "MiniBatchSize",64, ...
%     "InitialLearnRate",1e-3, ...
%     "ValidationData",{XVal2,YVal}, ...
%     "ValidationFrequency",25, ...
%     "ValidationPatience",10, ...
%     "Shuffle","every-epoch", ...
%     "ExecutionEnvironment", ternary(opts.useGPU,'gpu','cpu'), ...
%     "Verbose",opts.verbose, ...
%     "Plots","training-progress", ...
%     "CheckpointPath",opts.modelDir);
%
% net = trainNetwork(XTrain2,YTrain,lg,options);

%% ---------------------------- PLOTTING ---------------------------------
% if opts.savePlots
%     fig = figure('Name','Example epoch'); clf
%     plot(zscoref(XTrain{1})'); title('First normalised epoch'); xlabel('Time');
%     saveas(fig, fullfile(opts.modelDir,'example_epoch.png')); close(fig);
% end
% end % main

% ========================= HELPER FUNCTIONS ============================
function [cells,ref] = build_cells(files,doInit)
% Helper that converts each anot_data_all entry into a 2-D epoch matrix.
% Adds live feedback: files processed and cumulative epoch count.
if nargin<2, doInit = 0; end

cells      = cell(0);
ref        = [];
totalFiles = numel(files);
totEpochs  = 0;                       % running epoch counter

ft_progress('init','text', 'Reading %d files...', totalFiles);
for k = 1:totalFiles
    % -- load one .mat -------------------------------------------------
    S = load(files(k).name);
    if ~isfield(S,'anot_data_all'), continue; end
    
    % loop over trials in this file
    for t = 1:numel(S.anot_data_all)
        D = S.anot_data_all{t};
        if doInit && isempty(ref), ref = D.label; end
        D = do_ensure_consistent_labels(D, ref);
        allLbl  = D.label;
        megMask = startsWith(allLbl,'MEG');   % drop EEG, misc
        MEG_labels = allLbl(megMask);
        E  = D.trial{1}(megMask ,:);     % now 306 × T on Elekta
        E = maggrad_scale(E,MEG_labels);
        cells{end+1,1} = do_normalize_data(E,'demean'); %#ok<AGROW>
        totEpochs = totEpochs + 1;
    end
    % -- progress text ------------------------------------------------
    ft_progress(k/totalFiles, 'File %d/%d | epochs read: %d', ...
        k, totalFiles, totEpochs);
end
ft_progress('close');
fprintf('Finished build_cells: %d files, %d total epochs.',totalFiles, totEpochs);
end

% function out = parse_inputs(def,varargin)
% p = inputParser;
% addParameter(p,'ftRoot',def.ftRoot);
% addParameter(p,'codeRoot',def.codeRoot);
% addParameter(p,'spikeDir',def.spikeDir);
% addParameter(p,'noSpikeDir',def.noSpikeDir);
% addParameter(p,'skipList',def.skipList);
% addParameter(p,'useGPU',def.useGPU);
% addParameter(p,'modelDir',def.modelDir);
% addParameter(p,'savePlots',def.savePlots);
% addParameter(p,'verbose',def.verbose);
% addParameter(p,'valFrac',def.valFrac);
% parse(p,varargin{:});
% out = p.Results;
% end


function out = parse_inputs(def, varargin)
p = inputParser;
p.FunctionName = 'run_spike_classifier';
p.CaseSensitive = false;
p.PartialMatching = true;
p.KeepUnmatched = false;

% existing params
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

% NEW params for model switching
validArch = @(s) any(strcmpi(s,{'lstm','gru_td','transformer'}));
addParameter(p,'arch',def.arch, validArch);
addParameter(p,'dModel',def.dModel, @(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));
addParameter(p,'numHeads',def.numHeads, @(x) validateattributes(x,{'numeric'},{'scalar','integer','positive'}));

parse(p, varargin{:});
out = p.Results;

% extra sanity check
if strcmpi(out.arch,'transformer')
    assert(mod(out.dModel,out.numHeads)==0, ...
        'dModel (%d) must be divisible by numHeads (%d).', out.dModel, out.numHeads);
end
end


function setup_environment(o)
warning('off','all');
addpath(o.ftRoot); ft_defaults;
addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/functions_new')));
addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/helper')));
end

function r = ternary(c,a,b); if c, r=a; else, r=b; end; end


function X = maggrad_scale(X,chLabels)
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

