%% Spike-Detection MEG Toolkit (tidy, modular)
% Author: MCW MEG Lab  V. Youssofzadeh <vyoussofzadeh@mcw.edu>
% Updated: 2025-05-19
% -------------------------------------------------------------------------
% Top-level workflows:
%    run_spike_pipeline     PCA auto-encoder for QC / feature prep
%    run_spike_classifier   Bi-LSTM classifier (spike vs non-spike)
% Helper utilities are defined at the end of this file.
% -------------------------------------------------------------------------
% Example usage:
%   >> run_spike_pipeline
%   >> run_spike_classifier('useGPU',false,'verbose',true)
% -------------------------------------------------------------------------

%% -------------------------------------------------------------------------
function run_spike_DL_pipeline(varargin)
% Low-dimensional PCA trace + auto-regressive Bi-LSTM (unsupervised QC)

opts = default_opts();
opts = parse_inputs(opts,varargin{:});
setup_environment(opts);

files = rdir(fullfile(opts.dataDir,'annotated_data_anonymized','*.mat'));
assert(~isempty(files),'No annotated files found.');

[pcaRows, ~] = read_trials(files);            % cell array of [channels × time] matrices

% Pad / normalise
X = cellfun(@(x) pad_truncate(x,67), pcaRows,'uni',false);
[mu,sd] = get_zscore_stats(X);
X = cellfun(@(x) (x-mu)./sd,X,'uni',false);

% Model: auto-regressive forecast of next sample
layers = [
    sequenceInputLayer(67,'Normalization','zscore')
    bilstmLayer(128,'OutputMode','sequence')
    dropoutLayer(0.2)
    bilstmLayer(64,'OutputMode','last')
    fullyConnectedLayer(67)
    regressionLayer];

options = lstm_options(opts);
net = trainNetwork(X,X,layers,options);

save_model(net,options,'lstm_pca_qc',opts.modelRoot);
plot_mean_trace(X,opts,'mean_pca_qc');
end

%% -------------------------------------------------------------------------
function run_spike_classifier(varargin)
% Supervised Bi-LSTM spike detector (spike vs non-spike)

opts = default_opts();
opts.dataDirSpike   = fullfile(opts.dataDir,'annotated_data');
opts.dataDirNoSpike = fullfile(opts.dataDir,'annotated_data_nospike');
opts = parse_inputs(opts,varargin{:});
setup_environment(opts);

% ----- Load spike & non-spike epochs ------------------------------------
[spikeCells, refLabels]   = read_trials(rdir(fullfile(opts.dataDirSpike,'*.mat')));
[nospikeCells, ~]         = read_trials(rdir(fullfile(opts.dataDirNoSpike,'*.mat')),refLabels);

% Remove trials with NaNs
spikeCells   = spikeCells(~cellfun(@(x) any(isnan(x(:))),spikeCells));
nospikeCells = nospikeCells(~cellfun(@(x) any(isnan(x(:))),nospikeCells));

% Balance classes (optional)
minN = min(numel(spikeCells),numel(nospikeCells));
spikeCells   = spikeCells(1:minN);
nospikeCells = nospikeCells(1:minN);

% Build dataset & z-score each trial
X = [spikeCells; nospikeCells];
Y = categorical([ones(numel(spikeCells),1); zeros(numel(nospikeCells),1)], [0 1], {'NoSpike','Spike'});
X = cellfun(@(x) zscore(x,0,'all'), X,'uni',false);

% ----- Network ----------------------------------------------------------
inputSize = size(X{1},1);
numHidden1= 100; numHidden2 = 50;
layers = [
    sequenceInputLayer(inputSize)
    lstmLayer(numHidden1,'OutputMode','sequence')
    dropoutLayer(0.2)
    lstmLayer(numHidden2,'OutputMode','last')
    dropoutLayer(0.2)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

options = lstm_options(opts);
net = trainNetwork(X,Y,layers,options);

save_model(net,options,'lstm_spike_classifier',opts.modelRoot);
end

%% ------------------------- Helper Utilities ----------------------------
function opts = default_opts()
opts.ftRoot    = '/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
opts.codeRoot  = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git';
opts.dataDir   = '/data/MEG/Research/SpikeDectection/Epil_annotated_data';
opts.modelRoot = fullfile(opts.dataDir,'models');
opts.useGPU    = true;
opts.verbose   = true;
end

function opts = parse_inputs(opts,varargin)
% Merge user-supplied name/value pairs into opts struct
p = inputParser;
fn = fieldnames(opts);
for k = 1:numel(fn)
    addParameter(p,fn{k},opts.(fn{k}));
end
parse(p,varargin{:});
opts = p.Results;
end

function setup_environment(o)
warning('off','all'); restoredefaultpath;
addpath(o.ftRoot); ft_defaults;
addpath(genpath(fullfile(o.codeRoot,'FT_fucntions')));
end

function [cells, ref] = read_trials(fileList, ref)
if nargin<2, ref=[]; end
cells = {};
ft_progress('init','text','Reading trials...');
for k=1:numel(fileList)
    ft_progress(k/numel(fileList),'%d/%d',k,numel(fileList));
    S  = load(fileList(k).name,'anot_data_all');
    if isempty(ref); ref = S.anot_data_all{1}.label; end
    for t = 1:numel(S.anot_data_all)
        d = do_ensure_consistent_labels(S.anot_data_all{t},ref);
        x = do_normalize_data(d.trial{1}(:,1:67));
        cells{end+1,1} = pad_truncate(x,67); %#ok<AGROW>
    end
end
ft_progress('close');
end

function x = pad_truncate(x,N)
if size(x,2) < N
    x = [x, nan(size(x,1),N-size(x,2))];
else
    x = x(:,1:N);
end
end

function [mu,sd] = get_zscore_stats(X)
mu = mean(cat(2,X{:}),2);
sd = std(cat(2,X{:}),0,2);
end

function options = lstm_options(o)
options = trainingOptions('adam', ...
    'MaxEpochs',100, ...
    'MiniBatchSize',32, ...
    'InitialLearnRate',5e-3, ...
    'Shuffle','every-epoch', ...
    'Plots', ternary(o.verbose,'training-progress','none'), ...
    'Verbose', o.verbose, ...
    'ExecutionEnvironment', ternary(o.useGPU,'gpu','auto'), ...
    'CheckpointPath',o.modelRoot);
end

function save_model(net,opts,tag,root)
if ~exist(root,'dir'); mkdir(root); end
fname = fullfile(root,[tag,'_',datestr(now,'yyyymmdd_HHMM'),'.mat']);
save(fname,'net','opts'); fprintf('Model saved: %s\n',fname);
end

function plot_mean_trace(X,opts,tag)
if ~opts.verbose; return; end
fig = figure('Name',tag);
plot(mean(cat(2,X{:}),2)); xlabel('Time samples'); ylabel('Normalised PCA-1'); grid on;
if ~exist(opts.modelRoot,'dir'); mkdir(opts.modelRoot); end
saveas(fig,fullfile(opts.modelRoot,[tag,'.png']));
end

function r = ternary(c,a,b)
if c, r=a; else, r=b; end
end

%% End of file
