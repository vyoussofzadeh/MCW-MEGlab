%% Spike-Detection MEG Pipeline (tidy version)
% Deep-learning workflow for inter-ictal spike detection in MEG.
% Author: MCW MEG Lab  V. Youssofzadeh <vyoussofzadeh@mcw.edu>
% Last updated: 2025-05-19
% -------------------------------------------------------------------------
% Usage:
%   1. Edit the "USER SETTINGS" section below.
%   2. Run this script or call `run_spike_pipeline` from the command line.
% -------------------------------------------------------------------------
function run_spike_pipeline(varargin)

%% --------------------------- USER SETTINGS -----------------------------
opts.ftRoot     = '/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
opts.codeRoot   = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git';
% opts.dataDir    = '/data/MEG/Research/SpikeDectection/Epil_annotated_data/annotated_data_anonymized';
opts.dataDir    = '/data/MEG/Research/SpikeDectection/Epil_annotated_data/annotated_data';
opts.useGPU     = false;          % set false if GPU unavailable
opts.modelDir   = fullfile(opts.dataDir,'models');
opts.savePlots  = true;
opts.verbose    = true;

% Allow name-value pairs to override defaults
opts = parse_inputs(opts,varargin{:});

%% ----------------------------- ENV SETUP -------------------------------
setup_environment(opts);

%% -------------------------- LOAD ANNOTATED DATA ------------------------
files = rdir(fullfile(opts.dataDir,'*.mat'));
assert(~isempty(files),'No .mat files found in %s',opts.dataDir);

nFiles       = numel(files);
all_pca      = cell(nFiles,1);
refLabels    = [];

ft_progress('init','text','Loading and PCA-compressing trials...');
for f = 1:nFiles
    ft_progress(f/nFiles,'File %d of %d',f,nFiles);
    
    S = load(files(f).name);                % <-- loads everything
    if isfield(S, 'anot_data_all')
        trials = S.anot_data_all;
        % Initialise refLabels once
        if isempty(refLabels)
            refLabels = trials{1}.label;
        end
        % Concatenate PCA scores across trials for this run
        all_pca{f} = cellfun(@(d) pca_first_component(d,refLabels),trials,'UniformOutput',false);
    else
        [a,b] = fileparts(files(f).name);
        disp([b, ' was ignored, missing, annot_data_all'])
    end
end
ft_progress('close');

% Remove empty entries (files that were skipped)
all_pca = all_pca(~cellfun('isempty', all_pca));

% flatten nested cells: one entry per trial
X = [all_pca{:}];     % 1-by-N cell
X = X(:);             % column cell if you prefer

%%
%% ---------------------- PREPARE TRAINING DATA --------------------------
% Smooth and reshape every sequence to a 1×N row-vector
X = cellfun(@(x) reshape(smooth(x), 1, []), X, 'UniformOutput', false);

% Inspect lengths (optional)
lenVec = cellfun(@numel, X);
fprintf('Lengths: 67 -> %d   68 -> %d\n', sum(lenVec==67), sum(lenVec==68));

% ------- Global z-score normalisation --------
catAll = [X{:}];                 % 1 × totalTime
mu     = mean(catAll, 2);        % scalar
sd     = std( catAll, 0, 2 );
X      = cellfun(@(x) (x-mu)./sd, X, 'UniformOutput', false);
% ---------------------------------------------

%% --------------------------- LSTM ARCHITECTURE -------------------------
inputSize = 1;                   % row-vector per time step

layers = [
    sequenceInputLayer(inputSize,'Normalization','none')   % <-- fixed
    bilstmLayer(128,'OutputMode','sequence')
    dropoutLayer(0.2)
    bilstmLayer(64,'OutputMode','last')
    dropoutLayer(0.2)
    fullyConnectedLayer(inputSize)
    regressionLayer];

% Make sure the checkpoint folder exists
if ~exist(opts.modelDir, 'dir')
    mkdir(opts.modelDir);
end

options = trainingOptions('adam', ...
    'MaxEpochs',100, ...
    'MiniBatchSize',32, ...
    'InitialLearnRate',5e-3, ...
    'ValidationFrequency',25, ...
    'Shuffle','every-epoch', ...
    'Plots','training-progress', ...
    'Verbose',opts.verbose, ...
    'ExecutionEnvironment', ternary(opts.useGPU,'gpu','auto'), ...
    'CheckpointPath',opts.modelDir);

%% ----------------------------- TRAIN MODEL -----------------------------
Y   = X;                                     % auto-regressive target
net = trainNetwork(X, Y, layers, options);

if ~exist(opts.modelDir, 'dir'); mkdir(opts.modelDir); end
save(fullfile(opts.modelDir, 'lstm_spike_net.mat'), 'net', 'options', 'mu', 'sd');


%% ----------------------------- PLOTTING --------------------------------
if opts.savePlots
    fig = figure('Name','Mean PCA trace');
    plot(mean(cell2mat(X'),1)); xlabel('Time-points'); ylabel('Normalised PCA1'); grid on;
    saveas(fig,fullfile(opts.modelDir,'mean_pca.png'));
end

end % main function

% -------------------------- HELPER FUNCTIONS ----------------------------
function p = pca_first_component(anot_data,ref)
    % Ensure consistent channel order and return first PCA component (row-vector)
    anot_data = do_ensure_consistent_labels(anot_data,ref);
    comp      = do_pca(abs(anot_data.trial{1}),1);   % user-defined helper
    p         = comp(:)';                            % force row
end

function y = pad_truncate(x,N)
    % Pad with NaNs or truncate to length N
    if numel(x) < N, y = [x, nan(1,N-numel(x))]; else, y = x(1:N); end
end

function setup_environment(o)
    warning('off','all');
%     restoredefaultpath;
    addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/Squiggles/func')
    addpath(o.ftRoot); ft_defaults;
    addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/functions_new')));
    addpath(genpath(fullfile(o.codeRoot,'FT_fucntions/helper')));
end

function out = parse_inputs(def,varargin)
    p = inputParser;
    addParameter(p,'ftRoot',def.ftRoot);   %#ok<*DEP>
    addParameter(p,'codeRoot',def.codeRoot);
    addParameter(p,'dataDir',def.dataDir);
    addParameter(p,'useGPU',def.useGPU);
    addParameter(p,'modelDir',def.modelDir);
    addParameter(p,'savePlots',def.savePlots);
    addParameter(p,'verbose',def.verbose);
    parse(p,varargin{:});
    out = p.Results;
end

function r = ternary(cond,a,b)
    if cond, r=a; else, r=b; end
end
