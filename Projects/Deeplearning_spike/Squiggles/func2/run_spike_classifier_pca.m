%% Spike-Detection MEG Pipeline  Spike vs Non-Spike Classifier (July 22 2025)
% Deep-learning workflow to classify inter-ictal MEG epochs as **Spike** or **NoSpike**.
% Author: MCW MEG Lab  V. Youssofzadeh <vyoussofzadeh@mcw.edu>
% -------------------------------------------------------------------------
% Usage:
%   1. Edit the USER SETTINGS section.
%   2. Run the script or call `run_spike_classifier` from MATLAB.
% -------------------------------------------------------------------------
function run_spike_classifier_pca(varargin)

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
opts.valFrac     = 0.15;               % validation split
rng(0,'twister');                       % reproducible splits

% Allow name-value overrides
opts = parse_inputs(opts,varargin{:});

%% ----------------------------- ENV SETUP -------------------------------
setup_environment(opts);

gpuOK = false;
if opts.useGPU
    try, g=gpuDevice; gpuOK = g.ComputeCapability >= "3.0"; catch, end
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
[dataSpike, refLabels]   = proc_class(1, filesS);
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

%% ------------------------- NETWORK ARCHITECTURE ------------------------
inputSize = size(XTrain2{1},1);   % channels (rows)
numHidden1 = 128;  numHidden2 = 64;

layers = [
    sequenceInputLayer(inputSize,'Normalization','none')
    bilstmLayer(numHidden1,'OutputMode','sequence')
    dropoutLayer(0.2)
    bilstmLayer(numHidden2,'OutputMode','last')
    dropoutLayer(0.2)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',100, ...
    'MiniBatchSize',64, ...
    'InitialLearnRate',5e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',30, ...
    'LearnRateDropFactor',0.5, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{XVal2,YVal}, ...
    'ValidationFrequency',25, ...
    'ValidationPatience',10, ...
    'ExecutionEnvironment', ternary(opts.useGPU,'gpu','cpu'), ...
    'Verbose',opts.verbose, ...
    'Plots','training-progress', ...
    'CheckpointPath',opts.modelDir);

%% ------------------------------ TRAIN ----------------------------------
net = trainNetwork(XTrain2, YTrain, layers, options);

save(fullfile(opts.modelDir,'lstm_spike_classifier_pca.mat'), 'net','mu','sd','options');

%%
[YP,~] = classify(net,XVal2);
plotconfusion(YVal, YP);


%% ---------------------------- PLOTTING ---------------------------------
if opts.savePlots
    fig = figure('Name','Example epoch'); clf
    plot(zscoref(XTrain{1})'); title('First normalised epoch'); xlabel('Time');
    saveas(fig, fullfile(opts.modelDir,'example_epoch.png')); close(fig);
end
end % main



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
        E  = D.trial{1}(megMask ,1:67);     % now 306 × T on Elekta
        E = maggrad_scale(E,MEG_labels);
%         E =  do_pca(do_normalize_data(E,'demean'),1);
%         E1=[]; for j=1:size(E,1), E1(j,:) = smooth(smooth(E(j,:))); end
%         cells{end+1,1} = E1;
        cells{end+1,1} = E;
        totEpochs = totEpochs + 1;
    end
    
    X = cell2mat( cellfun(@(e) e(:).',  cells, 'UniformOutput',false) );  % N × (ch·T)

    % -- progress text ------------------------------------------------
    ft_progress(k/totalFiles, 'File %d/%d | epochs read: %d', ...
        k, totalFiles, totEpochs);
end
ft_progress('close');
fprintf('Finished build_cells: %d files, %d total epochs.',totalFiles, totEpochs);
end

function out = parse_inputs(def,varargin)
p = inputParser;
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
parse(p,varargin{:});
out = p.Results;
end

function setup_environment(o)
warning('off','all');
addpath(o.ftRoot); ft_defaults;
addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/functions_new')));
addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/helper')));
end

function r = ternary(c,a,b); if c, r=a; else, r=b; end; end

function unmixing = do_pca(sel_val, n)

numcomponent = n;
Nchans = size(sel_val,2);
dat = sel_val';
C = (dat*dat')./(size(dat,2)-1);

% eigenvalue decomposition (EVD)
[E,D] = eig(C);

% sort eigenvectors in descending order of eigenvalues
d = cat(2,(1:1:Nchans)',diag(D));
d = sortrows(d, -2);

% return the desired number of principal components
unmixing = E(:,d(1:numcomponent,1))';

end

function X = maggrad_scale(X,chLabels)
isMag  = endsWith(chLabels,'1');
isGrad = endsWith(chLabels,'2') | endsWith(chLabels,'3');
mMag  = median(abs(X(isMag ,:)), 'all','omitnan') + eps;
mGrad = median(abs(X(isGrad,:)), 'all','omitnan') + eps;
X(isMag ,:) = X(isMag ,:) ./ mMag;
X(isGrad,:) = X(isGrad,:) ./ mGrad;
end

%% ============================ HELPERS ==================================
%% Utility: stack epochs (channels × time × trials) and run global PCA
%   PCs = stack_pca(cells,k)
%    cells : cell array, each epoch ch×T (must be same size after pad)
%    k     : # of principal components to keep (default: all)
%   Returns PCs  (ch×T×k) so PC-i can be visualised as a channel×time map.
function PCs = stack_pca(cells,k)
    if nargin<2, k = min(5,numel(cells)); end
    ch = size(cells{1},1);  T = size(cells{1},2);
    cube = cat(3,cells{:});              % ch × T × N
    X    = reshape(cube, ch*T, [] ).';   % N × (ch·T)
    [coeff,~,~,~,~] = pca(X,'NumComponents',k);
    PCs = reshape(coeff(:,1:k), ch, T, k);   % ch × T × k
end
