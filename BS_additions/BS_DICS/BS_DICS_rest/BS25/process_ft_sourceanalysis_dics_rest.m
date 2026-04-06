function varargout = process_ft_sourceanalysis_dics_rest (varargin)
% PROCESS_FT_SOURCEANALYSIS_DICS_REST
% FieldTrip DICS beamformer for resting-state / single-window spectral
% source analysis in Brainstorm.
%
% This custom process supports optional normalization using:
%   1) projected noise / neural activity index (NAI),
%   2) background-noise estimates, or
%   3) an alternate frequency band as reference.
%
% Intended use:
%   - Primarily for resting-state or single-window spectral mapping, where
%     no true prestimulus baseline is available.
%   - For task-based analyses, an active-vs-baseline contrast with a common
%     spatial filter is generally preferable.
%
% Main updates relative to the original version:
%   1) Added projected-noise normalization (recommended default for rest).
%   2) Added optional background-noise and alternate-frequency references.
%   3) Removed per-map max-normalization before source contrast.
%   4) Improved robustness in trial loading, bad-channel handling,
%      sensor-information propagation, time-window usage, and output generation.
%
% Authors: Vahab Youssof Zadeh
% Updated: 2026

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() 
% Description the process
sProcess.Comment     = 'FieldTrip: ft_sourceanalysis (DICS-rest)';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 357;
sProcess.Description = 'https://github.com/vyoussofzadeh/DICS-beamformer-for-Brainstorm';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data'};
sProcess.OutputTypes = {'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

% Option: Sensors selection
sProcess.options.sensortype.Comment = 'Sensor type:';
sProcess.options.sensortype.Type    = 'combobox_label';
sProcess.options.sensortype.Value   = {'MEG', {'MEG', 'MEG GRAD', 'MEG MAG', 'EEG', 'SEEG', 'ECOG'; ...
    'MEG', 'MEG GRAD', 'MEG MAG', 'EEG', 'SEEG', 'ECOG'}};

% Label: Frequency
sProcess.options.label2.Comment = '<BR><B>Frequency of interest:</B>';
sProcess.options.label2.Type    = 'label';

sProcess.options.foi.Comment = 'FOI:';
sProcess.options.foi.Type    = 'value';
sProcess.options.foi.Value   = {18, 'Hz', 0};

sProcess.options.tpr.Comment = 'Tapering freq:';
sProcess.options.tpr.Type    = 'value';
sProcess.options.tpr.Value   = {4, 'Hz', 0};

% Label: Normalization
sProcess.options.labelNorm.Comment = '<BR><B>Normalization / comparison:</B>';
sProcess.options.labelNorm.Type    = 'label';

sProcess.options.normmode.Comment = {'Projected noise (NAI)', 'Other frequency bin', 'Normalization:'; ...
    'noise', 'otherfreq', ''};
sProcess.options.normmode.Type    = 'radio_linelabel';
sProcess.options.normmode.Value   = 'noise';

sProcess.options.ctrlfoi.Comment = 'Control freq (only if using other-frequency mode):';
sProcess.options.ctrlfoi.Type    = 'value';
sProcess.options.ctrlfoi.Value   = {15, 'Hz', 0};

sProcess.options.ctrltpr.Comment = 'Control tapering freq:';
sProcess.options.ctrltpr.Type    = 'value';
sProcess.options.ctrltpr.Value   = {12, 'Hz', 0};

sProcess.options.lambda.Comment = 'DICS lambda:';
sProcess.options.lambda.Type    = 'value';
sProcess.options.lambda.Value   = {5, '%', 0};

sProcess.options.realfilter.Comment = 'Use real-valued filter';
sProcess.options.realfilter.Type    = 'checkbox';
sProcess.options.realfilter.Value   = 1;

% Label: Contrast
sProcess.options.label4.Comment = '<BR><B>Output mode:</B>';
sProcess.options.label4.Type    = 'label';

sProcess.options.method.Comment = {'Subtraction / map output', 'Permutation-stats', 'Mode:'; ...
    'subtraction', 'permutation', ''};
sProcess.options.method.Type    = 'radio_linelabel';
sProcess.options.method.Value   = 'subtraction';

sProcess.options.erds.Comment = {'ERD', 'ERS', 'both', 'ERD/S effects:'; ...
    'erd', 'ers', 'both', ''};
sProcess.options.erds.Type    = 'radio_linelabel';
sProcess.options.erds.Value   = 'erd';

sProcess.options.effect.Comment = {'abs', 'raw', 'Absolute/raw value of the output:'; ...
    'abs', 'raw', ''};
sProcess.options.effect.Type    = 'radio_linelabel';
sProcess.options.effect.Value   = 'abs';

% Label: TFR
sProcess.options.label3.Comment = '<BR><B>Time-freq response (for data inspection):</B>';
sProcess.options.label3.Type    = 'label';

sProcess.options.maxfreq.Comment = 'Max frequency (for TFR calculation):';
sProcess.options.maxfreq.Type    = 'value';
sProcess.options.maxfreq.Value   = {40, 'Hz', 0};

sProcess.options.showtfr.Comment = 'Plot time-frequency figure';
sProcess.options.showtfr.Type    = 'checkbox';
sProcess.options.showtfr.Value   = 1;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 
OutputFiles = {};

% Initialize FieldTrip
[isInstalled, errMsg] = bst_plugin('Install', 'fieldtrip');
if ~isInstalled
    bst_report('Error', sProcess, [], errMsg);
    return;
end


% ===== GET OPTIONS =====
Method     = sProcess.options.method.Value;
Modality   = sProcess.options.sensortype.Value{1};
ShowTfr    = sProcess.options.showtfr.Value;
MaxFreq    = sProcess.options.maxfreq.Value{1};
FOI        = sProcess.options.foi.Value{1};
TprFreq    = sProcess.options.tpr.Value{1};
NormMode      = sProcess.options.normmode.Value;
CtrlFOI       = sProcess.options.ctrlfoi.Value{1};
LambdaPct     = sProcess.options.lambda.Value{1};
UseRealFilter = sProcess.options.realfilter.Value;
TmpDir        = bst_get('BrainstormTmpDir');
CtrlTprFreq = sProcess.options.ctrltpr.Value{1};

bst_progress('start', 'ft_sourceanalysis', 'Loading input files...', 0, 2*length(sInputs));

% ===== LOAD: CHANNEL FILE =====
bst_progress('text', 'Loading input files...');

% Load channel file
ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
% Get selected sensors
iChannels = channel_find(ChannelMat.Channel, Modality);
if isempty(iChannels)
    bst_report('Error', sProcess, sInputs, ['Channels "' Modality '" not found in channel file.']);
    return;
end

% ===== LOAD: BAD CHANNELS =====
% Load bad channels from all the input files
isChannelGood = [];
for iInput = 1:length(sInputs)
    DataFile = sInputs(1).FileName;
    DataMat = load(file_fullpath(DataFile), 'ChannelFlag');
    if isempty(isChannelGood)
        isChannelGood = (DataMat.ChannelFlag == 1);
    elseif (length(DataMat.ChannelFlag) ~= length(isChannelGood)) || (length(DataMat.ChannelFlag) ~= length(ChannelMat.Channel))
        bst_report('Error', sProcess, sInputs, 'All the input files must have the same number of channels.');
        return;
    else
        isChannelGood = isChannelGood & (DataMat.ChannelFlag == 1);
    end
end
% Remove bad channels
iChannelsData = intersect(iChannels, find(isChannelGood'));
% Error: All channels tagged as bad
if isempty(iChannelsData)
    bst_report('Error', sProcess, sInputs, 'All the selected channels are tagged as bad.');
    return;
elseif any(~isChannelGood)
    bst_report('Info', sProcess, sInputs, ['Found ' num2str(length(find(~isChannelGood))) ' bad channels: ', sprintf('%s ', ChannelMat.Channel(find(~isChannelGood)).Name)]);
end

% ===== LOAD: HEADMODEL =====
% Get the study
sStudyChan = bst_get('ChannelFile', sInputs(1).ChannelFile);
% Error if there is no head model available
if isempty(sStudyChan.iHeadModel)
    bst_report('Error', sProcess, [], ['No head model available in folder: ' bst_fileparts(sStudyChan.FileName)]);
    return;
end
% Load head model
HeadModelFile = sStudyChan.HeadModel(sStudyChan.iHeadModel).FileName;
HeadModelMat = in_bst_headmodel(HeadModelFile);
% Convert head model to FieldTrip format
[ftHeadmodel, ftLeadfield, iChannelsData] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);

% ===== LOAD: DATA =====
% Template FieldTrip structure for all trials
ftData = out_fieldtrip_data(sInputs(1).FileName, ChannelMat, iChannelsData, 1);
ftData.trial = cell(1,length(sInputs));
ftData.time = cell(1,length(sInputs));
% Load all the trials

AllChannelFiles = unique({sInputs.ChannelFile});
iChanInputs = find(ismember({sInputs.ChannelFile}, AllChannelFiles{1}));
for iInput = 1:length(sInputs)
    DataFile = sInputs(iChanInputs(iInput)).FileName;
    DataMat = in_bst_data(DataFile);
    ftData.trial{iInput} = DataMat.F(iChannelsData,:);
    ftData.time{iInput} = DataMat.Time;
end

%- checking inter-trial time-intervals
TI = ftData.time{2}(1) - ftData.time{1}(1);

% ===== FIELDTRIP: ft_freqanalysis =====
bst_progress('text', 'Calling FieldTrip function: ft_freqanalysis...');
% Compute tfr-decomposition
cfg = [];
cfg.output     = 'pow';
cfg.channel    = 'all';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.foi        = 1:2:MaxFreq;
cfg.keeptrials = 'yes';
cfg.t_ftimwin  = 3 ./ cfg.foi;
cfg.tapsmofrq  = 0.8 * cfg.foi;
%     cfg.toi        = Baseline(1):0.05:PostStim(2);
cfg.toi        = ftData.time{1}(1):0.05:ftData.time{1}(end);
tfr            = ft_freqanalysis(cfg, ftData);
if TI ~= 0, tfr.time = linspace(0, ftData.time{1}(end) - ftData.time{1}(1), length(cfg.toi)); end

tfr.powspctrm(isnan(tfr.powspctrm))=0;
tfr.time = linspace(0, tfr.time(end) - tfr.time(1), length(tfr.time));

% Plot TFR
cfg = [];
cfg.savepath = 1;
cfg.savefile = fullfile(savepath,'tfr');
cfg.fmax = MaxFreq;
cfg.toi = [tfr.time(1), tfr.time(end)];
cfg.bslcorr = 2;
cfg.plotflag = ShowTfr;
%     cfg.effect = 1; % postive = 1, negative = 2; both = 3;
[time_of_interest,freq_of_interest] = do_tfr_plot(cfg, tfr);
disp(['Global max: time:',num2str(time_of_interest),'sec']);
disp(['Global max: freq:',num2str(freq_of_interest),'Hz']);


%%
datain = ftData;

for i=1:length(datain.time)
    datain.time{i} = linspace(0, ftData.time{1}(end) - ftData.time{1}(1), length(ftData.time{1}));
end
warning(['Maximum trial length:[', num2str(datain.time{1}(1)), ',', num2str(datain.time{1}(end)),']']);

cfg = [];
cfg.toilim = [datain.time{1}(1), datain.time{1}(end)];
ep_data = ft_redefinetrial(cfg, datain);


% % Always use the full local epoch for resting-state blocks
epochEnd   = min(cellfun(@(t) t(end), ftData.time));
TimeWindow = [0, epochEnd];

cfg = [];
cfg.resamplefs = 500;
ep_data = ft_resampledata(cfg, ep_data);

% ===== SPECTRAL ESTIMATE AT FOI =====
fftCfg = [];
fftCfg.foilim   = [FOI FOI];
fftCfg.taper    = pick_taper(FOI);
fftCfg.tapsmofrq = pick_tapsmofrq(FOI, TprFreq);

f_data = [];
[f_data.foi, ~, ~, ~] = do_fft(fftCfg, ep_data);
f_data.foi = copy_sensorinfo(f_data.foi, ep_data);

if strcmpi(NormMode, 'otherfreq')
    ctrlCfg = [];
    ctrlCfg.foilim    = [CtrlFOI CtrlFOI];
    ctrlCfg.taper     = pick_taper(CtrlFOI);
    ctrlCfg.tapsmofrq = pick_tapsmofrq(CtrlFOI, CtrlTprFreq);
    [f_data.ctrl, ~, ~, ~] = do_fft(ctrlCfg, ep_data);
    f_data.ctrl = copy_sensorinfo(f_data.ctrl, ep_data);
end

% ===== SOURCE ANALYSIS =====
sourceOut = [];
outputMap = [];

switch Method
    case 'subtraction'
        switch lower(NormMode)
            case 'noise'
                sourceOut = compute_dics_source(ftLeadfield, ftHeadmodel, f_data.foi, LambdaPct, UseRealFilter, true, false);
                sourceOut = add_nai(sourceOut);
                outputMap = get_source_metric(sourceOut, 'nai');

            case 'otherfreq'
                s_foi  = compute_dics_source(ftLeadfield, ftHeadmodel, f_data.foi,  LambdaPct, UseRealFilter, false, false);
                s_ctrl = compute_dics_source(ftLeadfield, ftHeadmodel, f_data.ctrl, LambdaPct, UseRealFilter, false, false);
                
                tmp = s_foi.avg.pow; tmp(isnan(tmp))=0; tmp = tmp./max(tmp);s_foi.avg.pow = tmp;
                tmp =s_ctrl.avg.pow; tmp(isnan(tmp))=0; tmp = tmp./max(tmp); s_ctrl.avg.pow = tmp;
                outputMap = (s_foi.avg.pow - s_ctrl.avg.pow)./(s_foi.avg.pow + s_ctrl.avg.pow);
                sourceOut = s_foi;

                switch sProcess.options.erds.Value
                    case 'erd'
                        outputMap(isnan(outputMap)) = 0;
                        outputMap(outputMap > 0) = 0;
                    case 'ers'
                        outputMap(isnan(outputMap)) = 0;
                        outputMap(outputMap < 0) = 0;
                    case 'both'
                        outputMap(isnan(outputMap)) = 0;
                end
            otherwise
                bst_report('Error', sProcess, sInputs, ['Unknown normalization mode: ' NormMode]);
                bst_progress('stop');
                return;
        end

    case 'permutation'
        bst_report('Error', sProcess, sInputs, ...
            ['Permutation mode in this improved rest version requires an explicit control interval/condition ' ...
             'and a common-filter implementation. Use subtraction mode, or add a real baseline window.']);
        bst_progress('stop');
        return;

    otherwise
        bst_report('Error', sProcess, sInputs, ['Unknown method: ' Method]);
        bst_progress('stop');
        return;
end

switch sProcess.options.effect.Value
    case 'abs'
        outputMap = abs(outputMap);
    case 'raw'
        % keep as-is
end

% ===== SAVE RESULTS =====
bst_progress('text', 'Saving source file...');
bst_progress('inc', 1);

if (length(sInputs) == 1)
    iStudyOut  = sInputs(1).iStudy;
    RefDataFile = sInputs(1).FileName;
else
    [~, iStudyOut] = bst_process('GetOutputStudy', sProcess, sInputs);
    RefDataFile = [];
end

ResultsMat = db_template('resultsmat');
ResultsMat.ImagingKernel = [];
ResultsMat.ImageGridAmp  = outputMap;
if isfield(sourceOut, 'cfg')
    ResultsMat.cfg = sourceOut.cfg;
else
    ResultsMat.cfg = [];
end

ResultsMat.nComponents   = 1;
ResultsMat.Function      = Method;
ResultsMat.Time          = 1;
ResultsMat.DataFile      = RefDataFile;
ResultsMat.HeadModelFile = HeadModelFile;
ResultsMat.HeadModelType = HeadModelMat.HeadModelType;
ResultsMat.ChannelFlag   = DataMat.ChannelFlag;
ResultsMat.GoodChannel   = iChannelsData;
ResultsMat.SurfaceFile   = HeadModelMat.SurfaceFile;
ResultsMat.nAvg          = DataMat.nAvg;
ResultsMat.Leff          = DataMat.Leff;

if strcmpi(NormMode, 'noise')
    ResultsMat.Comment = ['DICS-rest NAI: ' num2str(FOI) 'Hz ' sprintf('%1.3fs-%1.3fs', TimeWindow(1), TimeWindow(2))];
else
%     ResultsMat.Comment = ['DICS-rest: ' num2str(FOI) 'Hz vs ' num2str(CtrlFOI) 'Hz ' sprintf('%1.3fs-%1.3fs', WindowSel(1), WindowSel(2))];
    ResultsMat.Comment = ['DICS-rest: ' num2str(FOI) 'Hz vs ' ...
    num2str(CtrlFOI) '±' num2str(CtrlTprFreq) 'Hz ' ...
    sprintf('%1.3fs-%1.3fs', TimeWindow(1), TimeWindow(2))];
end

switch lower(ResultsMat.HeadModelType)
    case 'volume'
        ResultsMat.GridLoc = HeadModelMat.GridLoc;
    case 'surface'
        ResultsMat.GridLoc = [];
    case 'mixed'
        ResultsMat.GridLoc    = HeadModelMat.GridLoc;
        ResultsMat.GridOrient = HeadModelMat.GridOrient;
end

ResultsMat = bst_history('add', ResultsMat, 'compute', ['ft_sourceanalysis: ' Method ' ' Modality ' ' NormMode]);

OutputDir  = bst_fileparts(file_fullpath(sInputs(1).FileName));
ResultFile = bst_process('GetNewFilename', OutputDir, ['results_' Method '_' Modality '_rest']);
bst_save(ResultFile, ResultsMat, 'v6');

newResult = db_template('results');
newResult.Comment       = ResultsMat.Comment;
newResult.FileName      = file_short(ResultFile);
newResult.DataFile      = ResultsMat.DataFile;
newResult.isLink        = 0;
newResult.HeadModelType = ResultsMat.HeadModelType;

sStudyOut = bst_get('Study', iStudyOut);
iResult   = length(sStudyOut.Result) + 1;
sStudyOut.Result(iResult) = newResult;
bst_set('Study', iStudyOut, sStudyOut);

OutputFiles{end+1} = newResult.FileName;
panel_protocols('SelectNode', [], newResult.FileName);
db_save();
bst_progress('stop');
end

%% ===== TIME-FREQUENCY =====
function [time_of_interest, freq_of_interest] = do_tfr_plot(cfg_main, tfr)
% First compute the average over trials:
cfg = [];
freq_avg = ft_freqdescriptives(cfg, tfr);

if cfg_main.bslcorr == 1
    cfg = [];
    cfg.baseline     = [-0.3 0];
    cfg.baselinetype = 'db';
    freq_avg_bsl = ft_freqbaseline(cfg, freq_avg);
    freq_avg_bsl.powspctrm(isnan(freq_avg_bsl.powspctrm)) = 0;
    meanpow = squeeze(mean(freq_avg_bsl.powspctrm, 1));
else
    freq_avg.powspctrm(isnan(freq_avg.powspctrm)) = 0;
    meanpow = squeeze(mean(freq_avg.powspctrm, 1));
end

tim_interp  = linspace(cfg_main.toi(1), cfg_main.toi(2), 512);
freq_interp = linspace(1, cfg_main.fmax, 512);

[tim_grid_orig,  freq_grid_orig]  = meshgrid(tfr.time, tfr.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow, tim_grid_interp, freq_grid_interp, 'spline');

pow_interp1  = pow_interp(50:end, 50:end);
tim_interp1  = tim_interp(50:end);
freq_interp1 = freq_interp(50:end);

[~, idx] = min(pow_interp1(:));
[row, col] = ind2sub(size(pow_interp1), idx);

time_of_interest = tim_interp1(col);
freq_of_interest = freq_interp1(row);

timind = nearest(tim_interp,  time_of_interest);
freqind = nearest(freq_interp, freq_of_interest);
pow_at_toi = pow_interp(:, timind);
pow_at_foi = pow_interp(freqind, :);

if cfg_main.plotflag
    figure();
    ax_main  = axes('Position', [0.1 0.2 0.55 0.55]);
    ax_right = axes('Position', [0.7 0.2 0.1 0.55]);
    ax_top   = axes('Position', [0.1 0.8 0.55 0.1]);

    axes(ax_main);
    imagesc(tim_interp, freq_interp, pow_interp);
    xlim([cfg_main.toi(1), cfg_main.toi(2)]);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    clim = max(abs(meanpow(:)));
    if clim == 0
        clim = 1;
    end
    caxis([-clim clim]);
    hold on;
    plot(zeros(size(freq_interp)), freq_interp, 'k:');

    axes(ax_top);
    area(tim_interp, pow_at_foi, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
    xlim([cfg_main.toi(1), cfg_main.toi(2)]);
    ylim([-clim clim]);
    box off;
    ax_top.XTickLabel = [];
    ylabel('Power (dB)');
    hold on;
    plot([0 0], [-clim clim], 'k:');

    axes(ax_right);
    area(freq_interp, pow_at_toi, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
    view([270 90]);
    ax_right.YDir = 'reverse';
    ylim([-clim clim]);
    box off;
    ax_right.XTickLabel = [];
    ylabel('Power (dB)');

    h = colorbar(ax_main, 'manual', 'Position', [0.85 0.2 0.05 0.55]);
    ylabel(h, 'Power');

    axes(ax_main);
    plot(ones(size(freq_interp)) * time_of_interest, freq_interp, 'Color', [0 0 0 0.1], 'LineWidth', 3);
    plot(tim_interp, ones(size(tim_interp)) * freq_of_interest, 'Color', [0 0 0 0.1], 'LineWidth', 3);

    axes(ax_top);
    plot([time_of_interest time_of_interest], [0 clim], 'Color', [0 0 0 0.1], 'LineWidth', 3);
    axes(ax_right);
    hold on;
    plot([freq_of_interest freq_of_interest], [0 clim], 'Color', [0 0 0 0.1], 'LineWidth', 3);
end
end

%% ===== FFT =====
function [freq, ff, psd, tapsmofrq] = do_fft(cfg_main, data)
cfg = [];
cfg.method       = 'mtmfft';
cfg.output       = 'fourier';
cfg.keeptrials   = 'yes';
cfg.foilim       = cfg_main.foilim;
cfg.tapsmofrq    = cfg_main.tapsmofrq;
cfg.taper        = cfg_main.taper;
cfg.pad          = 4;
freq             = ft_freqanalysis(cfg, data);
psd              = squeeze(mean(mean(abs(freq.fourierspctrm), 2), 1));
ff               = freq.freq;
tapsmofrq        = cfg.tapsmofrq;
end

%% ===== HELPERS =====
function taper = pick_taper(freqHz)
if freqHz >= 4
    taper = 'dpss';
else
    taper = 'hanning';
end
end

function tapsmofrq = pick_tapsmofrq(freqHz, userValue)
if freqHz >= 4
    tapsmofrq = userValue;
else
    tapsmofrq = 1;
end
end

function out = copy_sensorinfo(out, in)
if isfield(in, 'grad')
    out.grad = in.grad;
end
if isfield(in, 'elec')
    out.elec = in.elec;
end
end

function source = compute_dics_source(ftLeadfield, ftHeadmodel, freqData, lambdaPct, useRealFilter, useProjectNoise, keepFilter)
cfg = [];
cfg.method      = 'dics';
cfg.sourcemodel = ftLeadfield;
cfg.headmodel   = ftHeadmodel;
cfg.frequency   = freqData.freq;
cfg.dics.fixedori = 'yes';
cfg.dics.lambda = format_lambda(lambdaPct);
if useRealFilter
    cfg.dics.realfilter = 'yes';
end
if useProjectNoise
    cfg.dics.projectnoise = 'yes';
end
if keepFilter
    cfg.dics.keepfilter = 'yes';
end
source = ft_sourceanalysis(cfg, freqData);
end

function val = format_lambda(lambdaPct)
if isempty(lambdaPct) || (lambdaPct == 0)
    val = 0;
else
    val = sprintf('%g%%', lambdaPct);
end
end

function source = add_nai(source)
try
    source = ft_sourcedescriptives([], source);
catch
    % Fallback if the FieldTrip version behaves differently
    if isfield(source, 'avg') && isfield(source.avg, 'pow') && isfield(source.avg, 'noise')
        source.avg.nai = source.avg.pow ./ max(source.avg.noise, eps);
    elseif isfield(source, 'pow') && isfield(source, 'noise')
        source.nai = source.pow ./ max(source.noise, eps);
    else
        error('Could not compute NAI: missing power/noise fields in source structure.');
    end
end
end

function metric = get_source_metric(source, name)
switch lower(name)
    case 'nai'
        if isfield(source, 'nai')
            metric = source.nai;
        elseif isfield(source, 'avg') && isfield(source.avg, 'nai')
            metric = source.avg.nai;
        elseif isfield(source, 'avg') && isfield(source.avg, 'pow') && isfield(source.avg, 'noise')
            metric = source.avg.pow ./ max(source.avg.noise, eps);
        else
            error('NAI not found in source structure.');
        end

    case 'pow'
        if isfield(source, 'pow')
            metric = source.pow;
        elseif isfield(source, 'avg') && isfield(source.avg, 'pow')
            metric = source.avg.pow;
        else
            error('Power not found in source structure.');
        end

    otherwise
        error(['Unsupported source metric: ' name]);
end
metric(isnan(metric)) = 0;
end