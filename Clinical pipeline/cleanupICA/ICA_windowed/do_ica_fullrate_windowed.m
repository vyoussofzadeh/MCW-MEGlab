function [cleanData, badComponents, info] = do_ica_fullrate_windowed(cfg, data)
%DO_ICA_FULLRATE_WINDOWED Memory-bounded ICA cleanup without downsampling.
%
% ICA is learned from full-sampling-rate windows distributed across the
% recording. Only the selected component contributions are then subtracted
% from the complete recording, one block at a time. The original sampling
% rate and the signal outside the learned ICA subspace are preserved.
%
% Required cfg fields
%   layout             FieldTrip layout or layout filename
%
% Optional cfg fields
%   numcomponent       0/empty estimates rank (default 0)
%   nWindows           number of training windows (default 8)
%   windowSec          seconds per window (default 20)
%   maxAutoComponents  cap for automatically estimated rank (default 100)
%   rankTolerance      relative eigenvalue threshold (default 1e-8)
%   rankProbeSamples   samples used for rank estimate (default 50000)
%   cleanBlockSec      seconds cleaned per matrix block (default 20)
%   plotComponents     show component topographies (default true)
%   browseComponents   open ft_databrowser (default true)
%   badComponents      predefined components; otherwise prompt interactively

% This function requires one continuous FieldTrip trial.

% MCW pipeline helper, 2026.


arguments
    cfg struct
    data struct
end

cfg = set_default(cfg, 'numcomponent', 0);
cfg = set_default(cfg, 'nWindows', 8);
cfg = set_default(cfg, 'windowSec', 20);
cfg = set_default(cfg, 'maxAutoComponents', 100);
cfg = set_default(cfg, 'rankTolerance', 1e-8);
cfg = set_default(cfg, 'rankProbeSamples', 50000);
cfg = set_default(cfg, 'cleanBlockSec', 20);
cfg = set_default(cfg, 'plotComponents', true);
cfg = set_default(cfg, 'browseComponents', true);

if ~isfield(cfg, 'layout') || isempty(cfg.layout)
    error('do_ica_fullrate_windowed:MissingLayout', ...
        'cfg.layout is required.');
end
if ~isfield(data, 'trial') || numel(data.trial) ~= 1
    error('do_ica_fullrate_windowed:ContinuousDataRequired', ...
        'Expected one continuous FieldTrip trial.');
end
if isempty(data.trial{1}) || size(data.trial{1}, 2) < 2
    error('do_ica_fullrate_windowed:EmptyData', ...
        'The input recording contains no usable samples.');
end

fs = double(data.fsample);
nChannels = numel(data.label);
nSamples = size(data.trial{1}, 2);

fprintf('\nFull-rate windowed ICA\n');
fprintf('  Complete data: %d channels x %d samples at %.3f Hz\n', ...
    nChannels, nSamples, fs);

[rankEstimate, rankSpectrum] = estimate_rank(data.trial{1}, ...
    cfg.rankProbeSamples, cfg.rankTolerance);
autoComponents = min([rankEstimate, cfg.maxAutoComponents, nChannels]);

if isempty(cfg.numcomponent) || cfg.numcomponent <= 0
    nComponents = autoComponents;
    fprintf('  Estimated standardized rank: %d\n', rankEstimate);
    if autoComponents < rankEstimate
        fprintf('  Automatic component cap: %d\n', autoComponents);
    end
else
    nComponents = min(round(cfg.numcomponent), nChannels);
    fprintf('  User-requested ICA components: %d\n', nComponents);
end

if nComponents < 2
    error('do_ica_fullrate_windowed:RankTooLow', ...
        'At least two ICA components are required.');
end

[trainingData, trainingRanges] = make_training_data(data, ...
    cfg.nWindows, cfg.windowSec);
trainingSamples = sum(trainingRanges(:,2) - trainingRanges(:,1) + 1);
fprintf('  ICA training data: %d windows, %.1f seconds, %d samples\n', ...
    size(trainingRanges,1), trainingSamples/fs, trainingSamples);
fprintf('  Training retains the original %.3f-Hz sampling rate.\n', fs);
recommendedSamples = 30*nComponents^2;
if trainingSamples < recommendedSamples
    warning('do_ica_fullrate_windowed:LimitedTrainingSamples', ...
        ['ICA has %d training samples; approximately %d are recommended ' ...
         'for %d components. Increase nWindows/windowSec or reduce the ' ...
         'component count if convergence is unstable.'], ...
        trainingSamples, recommendedSamples, nComponents);
end

icaCfg = [];
icaCfg.method = 'runica';
icaCfg.numcomponent = nComponents;
icaCfg.feedback = 'text';

componentData = ft_componentanalysis(icaCfg, trainingData);

if cfg.plotComponents
    componentsPerFigure = 20;
    for firstComponent = 1:componentsPerFigure:nComponents
        lastComponent = min(firstComponent + componentsPerFigure - 1, ...
            nComponents);
        componentFigure = figure( ...
            'Name', sprintf('ICA components %d-%d', ...
                firstComponent, lastComponent), ...
            'NumberTitle', 'off', ...
            'Color', 'w');
        plotCfg = [];
        plotCfg.component = firstComponent:lastComponent;
        plotCfg.layout = cfg.layout;
        plotCfg.comment = 'no';
        figure(componentFigure);
        ft_topoplotIC(plotCfg, componentData);
        set(componentFigure, 'Name', sprintf('ICA components %d-%d', ...
            firstComponent, lastComponent), 'NumberTitle', 'off');
        drawnow;
    end
end

if cfg.browseComponents
    browseCfg = [];
    browseCfg.viewmode = 'component';
    browseCfg.layout = cfg.layout;
    browseCfg.blocksize = min(20, cfg.windowSec);
    browseCfg.continuous = 'no';
    ft_databrowser(browseCfg, componentData);
end

if isfield(cfg, 'badComponents') && ~isempty(cfg.badComponents)
    badComponents = cfg.badComponents;
else
    fprintf('\nReview the component figures/browser.\n');
    badComponents = input( ...
        'Enter artifact component numbers (e.g., [1 4 7], [] for none): ');
end

badComponents = unique(round(double(badComponents(:)')));
if any(~isfinite(badComponents)) || ...
        any(badComponents < 1 | badComponents > nComponents)
    error('do_ica_fullrate_windowed:InvalidComponents', ...
        'Artifact components must be integers from 1 through %d.', ...
        nComponents);
end

info = struct();
info.rankEstimate = rankEstimate;
info.rankSpectrum = rankSpectrum;
info.numcomponent = nComponents;
info.trainingRanges = trainingRanges;
info.trainingDurationSec = trainingSamples/fs;
info.fsample = fs;

if isempty(badComponents)
    fprintf('  No components selected; returning the original data.\n');
    cleanData = data;
    return;
end

% Copy only the small matrices needed for full-recording projection, then
% release the training data and component time courses before cleaning.
unmixing = double(componentData.unmixing(badComponents, :));
topography = double(componentData.topo(:, badComponents));
topoLabels = componentData.topolabel(:);

clear componentData trainingData

[labelsFound, dataChannelIndex] = ismember(topoLabels, data.label);
if ~all(labelsFound)
    missingLabels = strjoin(topoLabels(~labelsFound), ', ');
    error('do_ica_fullrate_windowed:ChannelMismatch', ...
        'ICA channels missing from complete data: %s', missingLabels);
end

% Assigning the function output back to the same caller variable allows
% recent MATLAB releases to reuse the input allocation where possible.
cleanData = data;
blockSamples = max(1, round(cfg.cleanBlockSec * fs));
nBlocks = ceil(nSamples / blockSamples);

fprintf('  Removing components %s from the complete recording...\n', ...
    mat2str(badComponents));
fprintf('  Cleaning in %d blocks of at most %.1f seconds.\n', ...
    nBlocks, blockSamples/fs);

for block = 1:nBlocks
    firstSample = (block-1)*blockSamples + 1;
    lastSample = min(block*blockSamples, nSamples);
    sampleIndex = firstSample:lastSample;

    originalClass = class(cleanData.trial{1});
    x = double(cleanData.trial{1}(dataChannelIndex, sampleIndex));
    blockMean = mean(x, 2);
    x = x - blockMean;

    artifactContribution = topography * (unmixing * x);
    x = x - artifactContribution + blockMean;

    cleanData.trial{1}(dataChannelIndex, sampleIndex) = ...
        cast(x, originalClass);

    if block == 1 || block == nBlocks || mod(block, 10) == 0
        fprintf('    block %d/%d\n', block, nBlocks);
    end
end

fprintf('  Full-rate component subtraction complete.\n');
end

function cfg = set_default(cfg, fieldName, defaultValue)
if ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
    cfg.(fieldName) = defaultValue;
end
end

function [rankEstimate, eigenvalues] = estimate_rank(dataMatrix, ...
        maximumSamples, relativeTolerance)
nSamples = size(dataMatrix, 2);
probeCount = min(nSamples, round(maximumSamples));
probeIndex = unique(round(linspace(1, nSamples, probeCount)));

x = double(dataMatrix(:, probeIndex));
x = x - mean(x, 2);
channelScale = sqrt(mean(x.^2, 2));
channelScale(~isfinite(channelScale) | channelScale <= 0) = 1;
x = x ./ channelScale;

correlationMatrix = (x*x') / max(1, size(x,2)-1);
correlationMatrix = (correlationMatrix + correlationMatrix') / 2;
eigenvalues = sort(real(eig(correlationMatrix)), 'descend');

if isempty(eigenvalues) || eigenvalues(1) <= 0
    rankEstimate = 0;
else
    rankEstimate = nnz(eigenvalues > eigenvalues(1)*relativeTolerance);
end
end

function [trainingData, ranges] = make_training_data(data, ...
        requestedWindows, requestedWindowSec)
fs = double(data.fsample);
nSamples = size(data.trial{1}, 2);
windowSamples = min(nSamples, max(2, round(requestedWindowSec*fs)));
maximumNonoverlappingWindows = max(1, floor(nSamples/windowSamples));
nWindows = min(max(1, round(requestedWindows)), ...
    maximumNonoverlappingWindows);

if nWindows == 1
    starts = max(1, floor((nSamples-windowSamples)/2) + 1);
else
    starts = round(linspace(1, nSamples-windowSamples+1, nWindows));
end

trainingData = data;
trainingData.trial = cell(nWindows, 1);
trainingData.time = cell(nWindows, 1);
trainingData.sampleinfo = zeros(nWindows, 2);
ranges = zeros(nWindows, 2);

for k = 1:nWindows
    firstSample = starts(k);
    lastSample = firstSample + windowSamples - 1;
    ranges(k,:) = [firstSample lastSample];

    trainingData.trial{k} = data.trial{1}(:, firstSample:lastSample);
    if isfield(data, 'time') && ~isempty(data.time)
        trainingData.time{k} = data.time{1}(firstSample:lastSample);
    else
        trainingData.time{k} = (0:windowSamples-1)/fs;
    end
    trainingData.sampleinfo(k,:) = [firstSample lastSample];
end
end
