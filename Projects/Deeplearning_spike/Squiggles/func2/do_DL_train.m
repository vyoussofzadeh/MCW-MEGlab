%% Spike-Detection MEG Pipeline (tidy version, REVISED 2025-07-22)
% Deep-learning workflow for inter-ictal spike detection in MEG.
% Author: MCW MEG Lab  V. Youssofzadeh <vyoussofzadeh@mcw.edu>
% -------------------------------------------------------------------------
% Usage:
%   1. Edit the "USER SETTINGS" section below.
%   2. Run this script or call `run_spike_pipeline` from the command line.
% -------------------------------------------------------------------------
function do_DL_train(varargin)

%% --------------------------- USER SETTINGS -----------------------------
opts.ftRoot     = '/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
opts.codeRoot   = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git';
opts.dataDir    = '/data/MEG/Research/SpikeDectection/Epil_annotated_data/annotated_data';
opts.modelDir   = fullfile(opts.dataDir,'models');
opts.useGPU     = false;          % auto-downgrades if incompatible
opts.savePlots  = true;
opts.verbose    = true;
opts.skipList   = {};            % e.g. {'pilot','bad_subject_12'}

% Allow name-value pairs to override defaults
opts = parse_inputs(opts,varargin{:});

%% ------------------------------ RNG -----------------------------------
rng(0,'twister');   % reproducible weights & data split

%% ----------------------------- ENV SETUP -------------------------------
setup_environment(opts);

%% --------------------------- LOAD DATA --------------------------------
files = rdir(fullfile(opts.dataDir,'*.mat'));
assert(~isempty(files),'No .mat files found in %s',opts.dataDir);

% --- optional skip list -------------------------------------------------
if ~isempty(opts.skipList)
    keep = true(size(files));
    for k = 1:numel(opts.skipList)
        keep = keep & ~contains({files.name}, opts.skipList{k});
    end
    files = files(keep);
end

nFiles    = numel(files);
all_pca   = cell(nFiles,1);
refLabels = [];

ft_progress('init','text','Loading & PCA-compressing trials...');
for f = 1:nFiles
    ft_progress(f/nFiles,'File %d of %d',f,nFiles);
    S = load(files(f).name);
    if isfield(S,'anot_data_all')
        trials = S.anot_data_all;
        if isempty(refLabels); refLabels = trials{1}.label; end
        all_pca{f} = cellfun(@(d) pca_first_component(d,refLabels), ...
                             trials,'UniformOutput',false);
    else
        warning('%s ignored  missing anot\_data\_all',files(f).name);
    end
end
ft_progress('close');

% Remove empties & flatten to one cell per trial
all_pca = all_pca(~cellfun('isempty',all_pca));
X = reshape([all_pca{:}],[],1);   % column cell

%% ---------------------- PREPARE TRAINING DATA --------------------------
X = cellfun(@(x) reshape(smooth(x),1,[]), X, 'UniformOutput',false);
lenVec = cellfun(@numel,X);
fprintf('Lengths: 67 -> %d   68 -> %d\n',sum(lenVec==67),sum(lenVec==68));

catAll = [X{:}];
mu = mean(catAll,2); sd = std(catAll,0,2);
X = cellfun(@(x)(x-mu)./sd,X,'UniformOutput',false);

%% -------------------------- DATA SPLIT ---------------------------------
cv = cvpartition(numel(X),'HoldOut',0.15);
idxTrain = training(cv);
idxVal   = test(cv);

XTrain = X(idxTrain);  YTrain = XTrain;   % seq-to-seq
XVal   = X(idxVal);    YVal   = XVal;

%% --------------------------- LSTM MODEL -------------------------------
inputSize = 1;
layers = [
    sequenceInputLayer(inputSize,'Normalization','none')
    bilstmLayer(128,'OutputMode','sequence')
    dropoutLayer(0.2)
    bilstmLayer(64,'OutputMode','sequence')
    dropoutLayer(0.2)
    fullyConnectedLayer(inputSize)
    regressionLayer];

% --- GPU fallback -------------------------------------------------------
if opts.useGPU
    try
        g = gpuDevice;
        if str2double(g.ComputeCapability) < 3.0
            warning('GPU compute capability <3.0  training on CPU.');
            opts.useGPU = false;
        end
    catch ME
        warning('No compatible GPU (%s)  training on CPU.',ME.message);
        opts.useGPU = false;
    end
end

exeEnv = ternary(opts.useGPU,'gpu','cpu');

% Ensure checkpoint directory exists
if ~exist(opts.modelDir,'dir'); mkdir(opts.modelDir); end

options = trainingOptions('adam', ...
    'MaxEpochs',100, ...
    'MiniBatchSize',64, ...
    'InitialLearnRate',5e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',30, ...
    'LearnRateDropFactor',0.5, ...
    'ValidationData',{XVal,YVal}, ...
    'ValidationFrequency',25, ...
    'ValidationPatience',10, ...
    'Shuffle','every-epoch', ...
    'Plots','training-progress', ...
    'Verbose',opts.verbose, ...
    'ExecutionEnvironment',exeEnv, ...
    'CheckpointPath',opts.modelDir);

%% ----------------------------- TRAIN -----------------------------------
net = trainNetwork(XTrain,YTrain,layers,options);
save(fullfile(opts.modelDir,'lstm_spike_net2.mat'),'net','options','mu','sd');

%% ----------------------------- PLOT ------------------------------------
if opts.savePlots
    fig = figure('Name','Mean PCA trace');
    plot(mean(cell2mat(X'),1)); xlabel('Time-points'); ylabel('Normalised PCA1'); grid on;
    saveas(fig,fullfile(opts.modelDir,'mean_pca.png'));
end

end % main

%% ----------------------- HELPER FUNCTIONS ------------------------------
function p = pca_first_component(anot_data,ref)
    anot_data = do_ensure_consistent_labels(anot_data,ref);
    comp      = do_pca(abs(anot_data.trial{1}),1);
    p         = comp(:)';
end

function setup_environment(o)
    warning('off','all');
    addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/Squiggles/func');
    addpath(o.ftRoot); ft_defaults;
    addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/functions_new')));
    addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/helper')));
end

function out = parse_inputs(def,varargin)
    p = inputParser;
    addParameter(p,'ftRoot',def.ftRoot);
    addParameter(p,'codeRoot',def.codeRoot);
    addParameter(p,'dataDir',def.dataDir);
    addParameter(p,'useGPU',def.useGPU);
    addParameter(p,'modelDir',def.modelDir);
    addParameter(p,'savePlots',def.savePlots);
    addParameter(p,'verbose',def.verbose);
    addParameter(p,'skipList',def.skipList);
    parse(p,varargin{:});
    out = p.Results;
end

function r = ternary(cond,a,b); if cond, r=a; else, r=b; end; end
