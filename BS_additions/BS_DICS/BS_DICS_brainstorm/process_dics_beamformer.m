function varargout = process_dics_beamformer(varargin)
% PROCESS_DICS_BEAMFORMER: Native Brainstorm DICS beamformer, no FieldTrip.
%
% The selected recordings are used to estimate a band-limited sensor CSD.
% Channels, active head model, leadfield, source grid, and optional noise
% covariance are obtained automatically from the Brainstorm database.
%
% Author:
%   Vahab YoussofZadeh, 2026

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() 
    sProcess.Comment     = 'DICS beamformer (Brainstorm)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 910;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/SourceEstimation';

    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'results', 'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.baseline.Comment = 'Baseline window:';
    sProcess.options.baseline.Type    = 'baseline';
    sProcess.options.baseline.Value   = [];

    sProcess.options.active.Comment = 'Active/post-stimulus window:';
    sProcess.options.active.Type    = 'poststim';
    sProcess.options.active.Value   = [];

    sProcess.options.freqrange.Comment = 'Frequency band:';
    sProcess.options.freqrange.Type    = 'range';
    sProcess.options.freqrange.Value   = {[8, 12], 'Hz', 2};

    sProcess.options.csd_win_length.Comment = 'CSD window length (0=longest shared length):';
    sProcess.options.csd_win_length.Type    = 'value';
    sProcess.options.csd_win_length.Value   = {0, 's', 3};

    sProcess.options.win_overlap.Comment = 'Window overlap:';
    sProcess.options.win_overlap.Type    = 'value';
    sProcess.options.win_overlap.Value   = {50, '%', 1};

    sProcess.options.sensortypes.Comment = 'Sensor types or names:';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG';

    sProcess.options.orientmode.Comment = { ...
        'Maximum-power orientation', ...
        'Head-model orientation', ...
        'NAI (uses study noise covariance)', ...
        'Source orientation:'; ...
        'max-power', 'fixed', 'nai', ''};
    sProcess.options.orientmode.Type  = 'radio_linelabel';
    sProcess.options.orientmode.Value = 'max-power';

    sProcess.options.contrast.Comment = { ...
        'Normalized difference: (post-pre)/(post+pre)', ...
        'Decibel ratio: 10*log10(post/pre)', ...
        'Raw difference: post-pre', ...
        'Post/pre ratio', ...
        'Output contrast:'; ...
        'normalized', 'db', 'difference', 'ratio', ''};
    sProcess.options.contrast.Type  = 'radio_linelabel';
    sProcess.options.contrast.Value = 'normalized';

    sProcess.options.reg.Comment = 'CSD regularization:';
    sProcess.options.reg.Type    = 'value';
    sProcess.options.reg.Value   = {0.05, 'fraction', 3};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    band = sProcess.options.freqrange.Value{1};
    Comment = sprintf('DICS contrast: %g-%g Hz', band(1), band(2));
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 
    OutputFiles = {};

    % Options
    baselineWindow = sProcess.options.baseline.Value{1};
    activeWindow = sProcess.options.active.Value{1};
    band = sProcess.options.freqrange.Value{1};
    requestedWinSec = sProcess.options.csd_win_length.Value{1};
    overlap = sProcess.options.win_overlap.Value{1};
    sensorTypes = strtrim(sProcess.options.sensortypes.Value);
    orientationMode = lower(sProcess.options.orientmode.Value);
    contrastMode = lower(sProcess.options.contrast.Value);
    regularization = sProcess.options.reg.Value{1};

    if numel(band) ~= 2 || any(~isfinite(band)) || band(1) < 0 || band(1) >= band(2)
        error('Frequency band must be [fmin fmax], with 0 <= fmin < fmax.');
    end
    if isempty(baselineWindow) || numel(baselineWindow) ~= 2 || baselineWindow(1) >= baselineWindow(2)
        error('Select a valid pre-stimulus baseline window.');
    end
    if isempty(activeWindow) || numel(activeWindow) ~= 2 || activeWindow(1) >= activeWindow(2)
        error('Select a valid active/post-stimulus window.');
    end
    if ~isscalar(requestedWinSec) || ~isfinite(requestedWinSec) || requestedWinSec < 0
        error('CSD window length must be zero (automatic) or a positive number of seconds.');
    end
    if ~isscalar(overlap) || ~isfinite(overlap) || overlap < 0 || overlap >= 100
        error('Window overlap must be between 0 and 100 percent.');
    end
    if ~isscalar(regularization) || ~isfinite(regularization) || regularization < 0
        error('Regularization must be a non-negative scalar.');
    end
    if isempty(sensorTypes)
        error('Select at least one sensor type.');
    end

    % All selected recordings must use the same channel file/head model.
    channelFile = sInputs(1).ChannelFile;
    if isempty(channelFile)
        error('The selected recording has no associated channel file.');
    end
    if any(~strcmp({sInputs.ChannelFile}, channelFile))
        error('All selected recordings must use the same channel file.');
    end
    ChannelMat = in_bst_channel(channelFile);
    nChannels = length(ChannelMat.Channel);

    % Exclude channels marked bad in any selected file.
    channelFlag = ones(nChannels, 1);
    for iFile = 1:length(sInputs)
        flagMat = in_bst_data(sInputs(iFile).FileName, 'ChannelFlag');
        if isfield(flagMat, 'ChannelFlag') && ~isempty(flagMat.ChannelFlag)
            if numel(flagMat.ChannelFlag) ~= nChannels
                error('Channel count mismatch in input file: %s', sInputs(iFile).FileName);
            end
            channelFlag(flagMat.ChannelFlag(:) < 0) = -1;
        end
    end
    goodChannels = good_channel(ChannelMat.Channel, channelFlag, sensorTypes);
    if isempty(goodChannels)
        error('No good channels match the sensor selection: %s.', sensorTypes);
    end

    % Locate the channel study and its active head model.
    [sStudyChannel, iStudyChannel] = bst_get('ChannelFile', channelFile);
    if isempty(sStudyChannel)
        iStudyChannel = sInputs(1).iStudy;
        sStudyChannel = bst_get('Study', iStudyChannel);
    end
    if isempty(sStudyChannel) || ~isfield(sStudyChannel, 'iHeadModel') || isempty(sStudyChannel.iHeadModel) || ...
            sStudyChannel.iHeadModel < 1 || sStudyChannel.iHeadModel > length(sStudyChannel.HeadModel)
        error('No active head model is available for the selected recordings.');
    end
    headModelFile = sStudyChannel.HeadModel(sStudyChannel.iHeadModel).FileName;
    HeadModel = in_bst_headmodel(headModelFile, 0, 'Gain', 'GridLoc', 'GridOrient', ...
        'GridAtlas', 'SurfaceFile', 'HeadModelType');
    if isempty(HeadModel.Gain) || size(HeadModel.Gain, 1) ~= nChannels
        error('The active head model does not match the selected channel file.');
    end
    if strcmpi(HeadModel.HeadModelType, 'mixed')
        error('Mixed head models are not supported by this DICS process.');
    end
    if mod(size(HeadModel.Gain, 2), 3) ~= 0
        error('Brainstorm head-model Gain must contain three columns per source.');
    end

    % Remove selected channels for which the forward model is invalid.
    validGainRows = all(isfinite(HeadModel.Gain), 2);
    goodChannels = goodChannels(validGainRows(goodChannels));
    if isempty(goodChannels)
        error('The head model has no valid gain rows for the selected sensors.');
    end

    % Apply active SSP projectors and the same EEG/ECOG/SEEG average-reference
    % transform to both recordings and leadfield.
    sensorTransform = eye(length(goodChannels));
    if isfield(ChannelMat, 'Projector') && ~isempty(ChannelMat.Projector)
        projector = process_ssp2('BuildProjector', ChannelMat.Projector, [1, 2]);
        if ~isempty(projector)
            sensorTransform = projector(goodChannels, goodChannels) * sensorTransform;
        end
    end
    selectedChannelTypes = {ChannelMat.Channel(goodChannels).Type};
    if any(ismember(unique(selectedChannelTypes), {'EEG', 'ECOG', 'SEEG'}))
        avgRef = panel_montage('GetMontageAvgRef', [], ChannelMat.Channel(goodChannels), ...
            channelFlag(goodChannels), 0);
        sensorTransform = avgRef.Matrix * sensorTransform;
    end

    Leadfield = sensorTransform * HeadModel.Gain(goodChannels, :);
    nSources = size(Leadfield, 2) / 3;

    % Estimate baseline and active CSDs separately. Both use the same actual
    % window length so their spectral smoothing and spatial filters are comparable.
    baselineCsdSum = [];
    activeCsdSum = [];
    nBaselineWindows = 0;
    nActiveWindows = 0;
    samplingRate = [];
    displayTime = [];
    effectiveWinSec = requestedWinSec;
    for iFile = 1:length(sInputs)
        bst_progress('text', sprintf('DICS: baseline/active CSD [%d/%d]...', iFile, length(sInputs)));
        [baselineData, baselineTime] = local_load_recording(sInputs(iFile), ChannelMat, goodChannels, baselineWindow);
        [activeData, activeTime] = local_load_recording(sInputs(iFile), ChannelMat, goodChannels, activeWindow);
        if isempty(baselineData) || size(baselineData, 2) < 2 || isempty(activeData) || size(activeData, 2) < 2
            error('No usable baseline or active samples were read from: %s', sInputs(iFile).FileName);
        end
        thisRate = 1 ./ median(diff(activeTime));
        baselineRate = 1 ./ median(diff(baselineTime));
        if abs(thisRate - baselineRate) > max(1e-6, thisRate * 1e-6)
            error('Baseline and active windows have inconsistent sampling rates.');
        end
        if isempty(samplingRate)
            samplingRate = thisRate;
            displayTime = [activeTime(1), activeTime(end)];
            if requestedWinSec == 0
                effectiveWinSec = min(size(baselineData, 2), size(activeData, 2)) ./ samplingRate;
            end
        elseif abs(thisRate - samplingRate) > max(1e-6, samplingRate * 1e-6)
            error('All selected recordings must have the same sampling rate.');
        end

        baselineData = sensorTransform * double(baselineData);
        activeData = sensorTransform * double(activeData);
        [baselineCsdSum, nAddedBaseline] = local_accumulate_csd( ...
            baselineCsdSum, baselineData, samplingRate, band, effectiveWinSec, overlap);
        [activeCsdSum, nAddedActive] = local_accumulate_csd( ...
            activeCsdSum, activeData, samplingRate, band, effectiveWinSec, overlap);
        nBaselineWindows = nBaselineWindows + nAddedBaseline;
        nActiveWindows = nActiveWindows + nAddedActive;
    end
    if nBaselineWindows < 1 || nActiveWindows < 1
        error('No baseline or active CSD windows were available. Check the selected windows.');
    end
    baselineCsd = baselineCsdSum ./ nBaselineWindows;
    activeCsd = activeCsdSum ./ nActiveWindows;
    baselineCsd = (baselineCsd + baselineCsd') ./ 2;
    activeCsd = (activeCsd + activeCsd') ./ 2;

    % One common filter is essential for an interpretable baseline/post contrast.
    commonCsd = (baselineCsdSum + activeCsdSum) ./ (nBaselineWindows + nActiveWindows);
    commonCsd = (commonCsd + commonCsd') ./ 2;

    solverOptions = {'Regularization', regularization, 'Orientation', orientationMode, ...
        'LeadfieldFormat', 'free'};

    if strcmp(orientationMode, 'fixed')
        fixedOrient = local_get_grid_orientations(HeadModel.GridOrient, nSources);
        solverOptions = [solverOptions, {'FixedOrient', fixedOrient}];
    elseif strcmp(orientationMode, 'nai')
        if ~isfield(sStudyChannel, 'NoiseCov') || isempty(sStudyChannel.NoiseCov)
            error('NAI requires a noise covariance in the study containing the channel file.');
        end
        noiseMat = load(file_fullpath(sStudyChannel.NoiseCov(1).FileName), 'NoiseCov');
        if ~isfield(noiseMat, 'NoiseCov') || ~isequal(size(noiseMat.NoiseCov), [nChannels, nChannels])
            error('The study noise covariance does not match the channel file.');
        end
        noiseCov = sensorTransform * noiseMat.NoiseCov(goodChannels, goodChannels) * sensorTransform';
        solverOptions = [solverOptions, {'NoiseCov', noiseCov}];
    end

    bst_progress('text', 'DICS: estimating source power...');
    [filters, ~, sourceOrientations, solverMeta] = ...
        bst_dics_beamformer(commonCsd, Leadfield, solverOptions{:});
    baselinePower = local_apply_filters(filters, baselineCsd);
    activePower = local_apply_filters(filters, activeCsd);
    sourceContrast = local_compute_contrast(activePower, baselinePower, contrastMode);

    % Create a standard Brainstorm source result. A static power map is
    % duplicated at the endpoints so Brainstorm can display it as a result.
    ResultsMat = db_template('resultsmat');
    ResultsMat.Comment = sprintf( ...
        'DICS %s %g-%g Hz | pre %g:%g ms | post %g:%g ms', ...
        contrastMode, band(1), band(2), ...
        1000 .* baselineWindow(1), 1000 .* baselineWindow(2), ...
        1000 .* activeWindow(1), 1000 .* activeWindow(2));
    ResultsMat.Function = 'DICS';
    ResultsMat.ImageGridAmp = [sourceContrast, sourceContrast];
    ResultsMat.ImagingKernel = [];
    ResultsMat.Time = displayTime;
    if length(sInputs) == 1
        ResultsMat.DataFile = sInputs(1).FileName;
    else
        ResultsMat.DataFile = [];
    end
    ResultsMat.HeadModelFile = headModelFile;
    ResultsMat.HeadModelType = HeadModel.HeadModelType;
    ResultsMat.ChannelFlag = channelFlag;
    ResultsMat.GoodChannel = goodChannels;
    ResultsMat.nComponents = 1;
    ResultsMat.nAvg = length(sInputs);
    if isfield(HeadModel, 'SurfaceFile')
        ResultsMat.SurfaceFile = file_short(HeadModel.SurfaceFile);
    end
    if strcmpi(HeadModel.HeadModelType, 'volume')
        ResultsMat.GridLoc = HeadModel.GridLoc;
    else
        ResultsMat.GridLoc = [];
    end
    ResultsMat.GridOrient = [];
    if isfield(HeadModel, 'GridAtlas')
        ResultsMat.GridAtlas = HeadModel.GridAtlas;
    end

    solverMeta.FrequencyBand = band;
    solverMeta.BaselineWindow = baselineWindow;
    solverMeta.ActiveWindow = activeWindow;
    solverMeta.BaselineWindowSec = baselineWindow;
    solverMeta.ActiveWindowSec = activeWindow;
    solverMeta.RequestedWindowLengthSec = requestedWinSec;
    solverMeta.EffectiveWindowLengthSec = effectiveWinSec;
    solverMeta.WindowOverlapPercent = overlap;
    solverMeta.WindowIntervalSec = effectiveWinSec .* (1 - overlap ./ 100);
    solverMeta.nBaselineCsdWindows = nBaselineWindows;
    solverMeta.nActiveCsdWindows = nActiveWindows;
    solverMeta.SamplingRate = samplingRate;
    solverMeta.SensorTypes = sensorTypes;
    solverMeta.ContrastMode = contrastMode;
    solverMeta.BaselinePower = baselinePower;
    solverMeta.ActivePower = activePower;
    solverMeta.SourceOrientations = sourceOrientations;
    ResultsMat.Options = solverMeta;
    ResultsMat = bst_history('add', ResultsMat, 'compute', ...
        sprintf(['DICS common-filter contrast: %g-%g Hz; baseline [%g, %g] s, ' ...
        'active [%g, %g] s; CSD window %g s, interval %g s (%g%% overlap); ' ...
        'baseline %d windows, active %d windows'], ...
        band(1), band(2), baselineWindow(1), baselineWindow(2), ...
        activeWindow(1), activeWindow(2), effectiveWinSec, ...
        solverMeta.WindowIntervalSec, overlap, nBaselineWindows, nActiveWindows));

    [sStudyOutput, iStudyOutput] = bst_process('GetOutputStudy', sProcess, sInputs);
    outputFile = bst_process('GetNewFilename', bst_fileparts(sStudyOutput.FileName), 'results_dics');
    bst_save(outputFile, ResultsMat, 'v6');
    db_add_data(iStudyOutput, outputFile, ResultsMat);
    OutputFiles = {outputFile};
end


%% ===== APPLY COMMON FILTERS =====
function sourcePower = local_apply_filters(filters, Csd)
    sourcePower = real(sum(conj(filters) .* (Csd * filters), 1))';
    finiteValues = sourcePower(isfinite(sourcePower));
    if isempty(finiteValues)
        return;
    end
    powerScale = max(abs(finiteValues));
    if powerScale > 0
        numericalTolerance = 10 .* eps(class(sourcePower)) .* powerScale;
    else
        numericalTolerance = 0;
    end
    sourcePower(sourcePower < -numericalTolerance) = NaN;
    sourcePower(sourcePower < 0 & sourcePower >= -numericalTolerance) = 0;
end


%% ===== COMPUTE BASELINE/ACTIVE CONTRAST =====
function sourceContrast = local_compute_contrast(activePower, baselinePower, contrastMode)
    finitePower = [activePower(isfinite(activePower)); baselinePower(isfinite(baselinePower))];
    positivePower = finitePower(finitePower > 0);
    if isempty(positivePower)
        sourceContrast = NaN(size(activePower));
        return;
    end
    powerScale = max(positivePower);
    powerFloor = max(realmin(class(activePower)), eps(class(activePower)) .* powerScale);
    switch contrastMode
        case 'normalized'
            denominator = activePower + baselinePower;
            sourceContrast = (activePower - baselinePower) ./ denominator;
            invalid = ~isfinite(activePower) | ~isfinite(baselinePower) | denominator <= powerFloor;
            sourceContrast(invalid) = NaN;
        case 'db'
            activeSafe = max(activePower, powerFloor);
            baselineSafe = max(baselinePower, powerFloor);
            sourceContrast = 10 .* log10(activeSafe ./ baselineSafe);
            sourceContrast(~isfinite(activePower) | ~isfinite(baselinePower)) = NaN;
        case 'difference'
            sourceContrast = activePower - baselinePower;
        case 'ratio'
            activeSafe = max(activePower, powerFloor);
            baselineSafe = max(baselinePower, powerFloor);
            sourceContrast = activeSafe ./ baselineSafe;
            sourceContrast(~isfinite(activePower) | ~isfinite(baselinePower)) = NaN;
        otherwise
            error('Unknown DICS contrast mode: %s.', contrastMode);
    end
end


%% ===== LOAD ONE RECORDING =====
function [F, Time] = local_load_recording(sInput, ChannelMat, goodChannels, timeWindow)
    if strcmpi(sInput.FileType, 'raw')
        DataMat = in_bst_data(sInput.FileName, 'F', 'Time');
        sFile = DataMat.F;
        if isempty(timeWindow)
            sampleBounds = [];
        else
            timeIndices = bst_closest(timeWindow, DataMat.Time);
            sampleBounds = round(sFile.prop.times(1) .* sFile.prop.sfreq) + timeIndices - 1;
        end
        [F, Time] = in_fread(sFile, ChannelMat, 1, sampleBounds, goodChannels);
    else
        DataMat = in_bst_data(sInput.FileName, 'F', 'Time');
        if ~isnumeric(DataMat.F)
            error('Input is not a numeric Brainstorm data matrix: %s', sInput.FileName);
        end
        F = DataMat.F(goodChannels, :);
        Time = DataMat.Time;
        if ~isempty(timeWindow)
            timeIndices = bst_closest(timeWindow, Time);
            sampleIndices = timeIndices(1):timeIndices(end);
            F = F(:, sampleIndices);
            Time = Time(sampleIndices);
        end
    end
    Time = double(Time(:)');
end


%% ===== ACCUMULATE CSD =====
function [CsdSum, nWindows] = local_accumulate_csd(CsdSum, F, samplingRate, band, winSec, overlap)
    if band(2) > samplingRate ./ 2
        error('Upper frequency (%g Hz) exceeds Nyquist (%g Hz).', band(2), samplingRate ./ 2);
    end

    [nChannels, nSamples] = size(F);
    winSamples = max(2, round(winSec .* samplingRate));
    if winSamples > nSamples
        error('CSD window (%g s) is longer than the selected data segment (%g s).', ...
            winSec, nSamples ./ samplingRate);
    end
    step = max(1, round(winSamples .* (1 - overlap ./ 100)));
    nfft = 2 ^ nextpow2(winSamples);
    frequencies = (0:floor(nfft ./ 2)) .* samplingRate ./ nfft;
    bandIndices = find(frequencies >= band(1) & frequencies <= band(2));
    if isempty(bandIndices)
        error('No FFT bins fall inside %g-%g Hz. Increase the CSD window length.', band(1), band(2));
    end

    % Periodic Hann window, implemented locally to avoid toolbox dependencies.
    window = 0.5 - 0.5 .* cos(2 .* pi .* (0:winSamples-1)' ./ winSamples);
    windowPower = sum(window .^ 2);
    oneSidedWeight = 2 .* ones(1, length(bandIndices));
    oneSidedWeight(frequencies(bandIndices) == 0) = 1;
    if rem(nfft, 2) == 0
        oneSidedWeight(frequencies(bandIndices) == samplingRate ./ 2) = 1;
    end

    if isempty(CsdSum)
        CsdSum = complex(zeros(nChannels, nChannels));
    elseif ~isequal(size(CsdSum), [nChannels, nChannels])
        error('Selected recordings produced inconsistent channel counts.');
    end

    nWindows = 0;
    for iStart = 1:step:(nSamples - winSamples + 1)
        segment = F(:, iStart:(iStart + winSamples - 1));
        segment = bsxfun(@minus, segment, mean(segment, 2));
        segment = bsxfun(@times, segment, window');
        spectrum = fft(segment, nfft, 2);
        spectrum = spectrum(:, bandIndices);
        spectrum = bsxfun(@times, spectrum, sqrt(oneSidedWeight));
        CsdSum = CsdSum + (spectrum * spectrum') ./ ...
            (samplingRate .* windowPower .* length(bandIndices));
        nWindows = nWindows + 1;
    end
end


%% ===== GRID ORIENTATIONS =====
function fixedOrient = local_get_grid_orientations(gridOrient, nSources)
    if isempty(gridOrient)
        error('Head-model orientation was selected, but GridOrient is empty.');
    end
    if isequal(size(gridOrient), [nSources, 3])
        fixedOrient = gridOrient';
    elseif isequal(size(gridOrient), [3, nSources])
        fixedOrient = gridOrient;
    else
        error('GridOrient dimensions do not match the number of sources.');
    end
end
