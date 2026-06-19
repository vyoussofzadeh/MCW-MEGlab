function Step7a_Interactive_Review_Completed_Model_Events()
%% ------------------------------------------------------------------------
% Step7_Interactive_Review_Completed_Model_Events.m
%
% Interactive review/feedback of completed model-detected events.
%
% Input:
%   Step 6 full MAT files:
%     *_model_candidates_gt75_refined_full.mat
%
% For each event, displays:
%   1) Butterfly plot around event
%   2) GFP/global energy trace
%   3) Sensor-layout waveform map
%   4) Spatial sensor amplitude map at local GFP peak
%
% Keyboard:
%   1 = true spike
%   0 = false positive / no-spike
%   u = unsure
%   n/right/space = next without changing label
%   b/left = back
%   s = save
%   q/escape = save and quit
%
% Output per Step 6 file:
%   *_interactive_reviewed.mat
%   *_interactive_reviewed.csv
%   *_interactive_reviewed.xlsx
%
% Output table has:
%   reviewLabel:
%       1   true spike
%       0   false positive / no-spike
%       NaN unsure / not reviewed
%
%   reviewType:
%       true_spike
%       false_positive
%       unsure
% -------------------------------------------------------------------------

clc;

%% ======================== USER SETTINGS =================================

opts.ftRoot  = '/home/vyoussofzadeh/Desktop/research_workspace/tools/fieldtrip/fieldtrip_2022/';
opts.mneRoot = '/MEG_data/MEG_Tools/mne';

% % Folder containing Step 6 completed full MAT files
% opts.step6OutDir = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data/model_outputs_gt75_bp5_50_m250_p500';
% % Pattern from Step 6 batch output
% opts.filePattern = '*_model_candidates_gt75_refined_full.mat';

opts.step6OutDir = '/home/vyoussofzadeh/github/MCW-MEGlab/MCW-MEGlab/Projects/Deeplearning_spike/Data/';
opts.filePattern = 'baer_jessica_Run02_spont_eyesclosed_raw_t_sss_ecgClean_raw_DS_model_candidates_gt70_refined_full.mat';

% opts.navigatorTableFontSize = 5;   % try 6 if still too large

% Review window around each event
opts.reviewWinSec = [-0.750 0.750];

% Local peak search window shown in display
opts.peakSearchSec = [-0.150 0.150];

% Display options
opts.nButterflyChannels = 40;     % top active channels in butterfly
opts.maxSensorWaveforms = 306;    % lower to 120 if display is slow

opts.sensorMapSensorType = 'mag';

% Apply same filter during display if Step 6 opts says so
opts.matchStep6Bandpass = true;

% Continue where you left off if reviewed file exists
opts.resumeExistingReview = true;

% Export reviewed EVE file after save
opts.exportReviewedEve = true;
opts.reviewedSpikeCode = 1111;
opts.reviewedNoSpikeCode = 2222;

% Review navigation
opts.showReviewNavigator = true;
% opts.navigatorTablePosition = [0.03 0.03 0.94 0.91];
opts.navigatorTablePosition = [0.02 0.03 0.55 0.92];
% opts.navigatorColumnWidth = {35, 65, 65, 80, 95};
opts.navigatorColumnWidth = {55, 95, 95, 110, 170};

% Thin probability strip
opts.axProbPos = [0.61 0.875 0.34 0.045];

% Move GFP lower and make it a bit smaller
opts.axGfpPos  = [0.61 0.620 0.34 0.120];

% Move sensor map lower if needed, but keep it large
opts.axMapPos  = [0.55 0.035 0.42 0.530];

% Candidate table uses the full right panel
opts.navigatorTablePosition = [0.03 0.03 0.94 0.91];

opts.showProbabilityNavigator = true;
opts.probNavigatorShowThreshold = true;
opts.probNavigatorBarWidth = 1.0;

% Probability bar plot display
opts.probNavigatorBarWidth = 0.45;   % smaller = more space between bars
opts.probNavigatorXLabelEvery = 1;   % show every candidate number
opts.probNavigatorEdgeColor = [0.15 0.15 0.15];
opts.probNavigatorLineWidth = 0.25;

%% ---- Probability navigator ---------------------------------------------

opts.showProbabilityNavigator = true;

% Right-side navigator layout
opts.probNavigatorAxesPosition = [0.06 0.67 0.88 0.27];  % top bar plot
opts.navigatorTablePosition    = [0.03 0.03 0.94 0.58];  % bottom table

% Probability display
opts.probNavigatorShowThreshold = true;
opts.probNavigatorBarWidth = 1.0;


%% ---- Single-window GUI layout ------------------------------------------

opts.mainFigPosition    = [0.1 0.2 0.8 0.6];

% Top buttons
opts.controlPanelPos    = [0.00 0.89 1.00 0.11];

% Full-width probability strip below buttons
opts.probPanelPos       = [0.00 0.82 1.00 0.07];

% Main lower panels
opts.plotPanelPos       = [0.00 0.00 0.70 0.82];
opts.navigatorPanelPos  = [0.70 0.00 0.30 0.82];

% Axes inside left plot panel only
opts.axTimePos = [0.08 0.08 0.48 0.86];
opts.axGfpPos  = [0.62 0.73 0.34 0.15];
opts.axMapPos  = [0.56 0.05 0.41 0.60];

% Candidate table fills right panel
opts.navigatorTablePosition = [0.02 0.03 0.96 0.92];
opts.navigatorTableFontSize = 9;
% opts.navigatorColumnWidth = {35, 75, 75, 105, 110};
% opts.navigatorColumnWidth = {28, 60, 55, 75, 80};

% Probability strip display
opts.probNavigatorBarWidth = 0.35;
opts.probNavigatorXLabelEvery = 1;     % show 1:n candidates
opts.probNavigatorXTickRotation = 90;
opts.probNavigatorEdgeColor = [0.10 0.10 0.10];
opts.probNavigatorLineWidth = 0.20;
opts.probNavigatorShowThreshold = true;

%% ---- Improved review display options -----------------------------------

% Load ECG along with MEG for visual inspection
opts.reviewLoadChannels = {'MEG','ECG'};

% Main left time-course panel
% Options: 'mag', 'grad', 'meg', 'mag_ecg', 'meg_ecg', 'ecg', 'all'
opts.timeCourseSensorType = 'mag_ecg';

% Time-course scaling options
% 'per_channel' is best for inspection because MAG and ECG have different units.
% Options: 'per_channel', 'global', 'manual'
opts.timeCourseScaleMode = 'per_channel';
opts.timeCourseScalePrctile = 98;
opts.timeCourseManualScale = [];     % used only if mode = 'manual'
opts.timeCourseAmpMultiplier = 1.0;
opts.timeCourseOffset = 2.5;

% How many traces to show in the left panel
% For all magnetometers + ECG, this should be about 102 + ECG.
opts.maxTimeCourseChannels = 130;

% Sensor waveform map: show only this interval around spike time
opts.sensorMapWinSec = [-0.150 0.150];

% Sensor waveform map mini-trace scaling
opts.maxSensorWaveforms = 306;
opts.sensorMapTimeScale = 0.12;
opts.sensorMapAmpScale  = 0.045;

% Neuromag layout
opts.layoutFile = fullfile(opts.ftRoot, 'template', 'layout', 'neuromag306all.lay');
if ~isfile(opts.layoutFile)
    opts.layoutFile = 'neuromag306all.lay';
end

% Panel layout positions: [left bottom width height]

% Left time-course panel
opts.axTimePos = [0.085 0.08 0.45 0.82];

% Make GFP shorter/smaller
opts.axGfpPos  = [0.61 0.76 0.34 0.10];

% Make sensor waveform/layout panel bigger
opts.axMapPos  = [0.54 0.05 0.44 0.64];

% Timecourse y-axis sensor labels
opts.timeCourseShowYLabels = true;
opts.timeCourseYLabelEvery = 1;      % label every 5th sensor to avoid crowding
opts.timeCourseLabelFontSize = 6;

% Sensor waveform map labels
opts.sensorMapShowLabels = true;
opts.sensorMapLabelEvery = 12;        % use 1 for every sensor, but may be crowded
opts.sensorMapLabelFontSize = 5;

%% ---- Improved display layout settings ----------------------------------

% Use Neuromag/Elekta 306 layout if available
opts.layoutFile = fullfile(opts.ftRoot, 'template', 'layout', 'neuromag306all.lay');

% If the full path does not exist, FieldTrip may still find it by name
if ~isfile(opts.layoutFile)
    opts.layoutFile = 'neuromag306all.lay';
end

% Large left butterfly panel
opts.nButterflyChannels = 80;

% Right-bottom mbrowse-like sensor waveform map
opts.maxSensorWaveforms = 120;
opts.sensorMapTimeScale = 0.08;   % horizontal width of mini waveforms
opts.sensorMapAmpScale  = 0.030;  % vertical height of mini waveforms

%% ======================== INITIALIZE ====================================

if exist(opts.ftRoot, 'dir')
    addpath(opts.ftRoot);
    ft_defaults;
else
    error('FieldTrip folder not found: %s', opts.ftRoot);
end

if exist(opts.mneRoot, 'dir')
    addpath(genpath(opts.mneRoot));
end

rehash toolboxcache;

%% ======================== FIND COMPLETED EVENT FILES =====================

L = dir(fullfile(opts.step6OutDir, opts.filePattern));

if isempty(L)
    error('No Step 6 full MAT files found:\n%s', fullfile(opts.step6OutDir, opts.filePattern));
end

step6Files = fullfile({L.folder}, {L.name})';

fprintf('\nFound %d completed Step 6 files:\n', numel(step6Files));
for i = 1:numel(step6Files)
    fprintf('[%03d] %s\n', i, step6Files{i});
end

sel = 1:numel(step6Files);

%% ======================== REVIEW SELECTED FILES ==========================

for ff = sel(:)'

    fprintf('\n============================================================\n');
    fprintf('Reviewing file %d/%d:\n%s\n', ff, numel(step6Files), step6Files{ff});
    fprintf('============================================================\n');

    review_one_step6_file(step6Files{ff}, opts);
end

fprintf('\nAll selected files reviewed/saved.\n');

end

%% ========================================================================
% Main review function
% ========================================================================

function review_one_step6_file(step6MatFile, opts)

[step6Dir, step6Base, ~] = fileparts(step6MatFile);

outMat  = fullfile(step6Dir, [step6Base '_interactive_reviewed.mat']);
outCsv  = fullfile(step6Dir, [step6Base '_interactive_reviewed.csv']);
outXlsx = fullfile(step6Dir, [step6Base '_interactive_reviewed.xlsx']);

%% ---- Load Step 6 output ----

S = load(step6MatFile);

if ~isfield(S, 'Tfull')
    error('Step 6 MAT file does not contain Tfull: %s', step6MatFile);
end

if ~isfield(S, 'opts')
    error('Step 6 MAT file does not contain opts: %s', step6MatFile);
end

Tfull = S.Tfull;
opts6 = S.opts;

if height(Tfull) == 0
    warning('No events in Tfull: %s', step6MatFile);
    return;
end

rawFile = opts6.rawFile;

if ~isfile(rawFile)
    error('Raw FIF file not found:\n%s', rawFile);
end

% Resume previous labels if available
if opts.resumeExistingReview && isfile(outMat)
    R = load(outMat, 'Treviewed');
    if isfield(R, 'Treviewed') && height(R.Treviewed) == height(Tfull)
        fprintf('Resuming existing review:\n%s\n', outMat);
        Tfull = R.Treviewed;
    end
end

% Ensure review columns exist
if ~ismember('reviewLabel', Tfull.Properties.VariableNames)
    Tfull.reviewLabel = nan(height(Tfull),1);
end

if ~ismember('reviewType', Tfull.Properties.VariableNames)
    Tfull.reviewType = repmat({''}, height(Tfull), 1);
end

if ~ismember('reviewNotes', Tfull.Properties.VariableNames)
    Tfull.reviewNotes = repmat({''}, height(Tfull), 1);
end

%% ---- Load MEG ----

fprintf('\nLoading raw MEG:\n%s\n', rawFile);

[dataMat, chanLabels, Fs] = load_meg_for_review(rawFile, opts, opts6);

% Identify channel groups for display
chInfo = get_channel_groups(chanLabels);

if numel(chInfo.megIdx) < 306
    warning('Fewer than 306 MEG channels found: %d', numel(chInfo.megIdx));
end

% Model/GFP/sensor map should use MEG channels only
nMegUse = min(306, numel(chInfo.megIdx));
megModelIdx = chInfo.megIdx(1:nMegUse);

% If model normalization has Cmodel channels, use exactly that many MEG channels
if exist('Cmodel','var') && numel(megModelIdx) >= Cmodel
    megModelIdx = megModelIdx(1:Cmodel);
end

fprintf('Loaded MEG: %d channels x %d samples, Fs=%.3f Hz\n', ...
    size(dataMat,1), size(dataMat,2), Fs);

% Use model Fs if available
if isfield(S, 'FsModel')
    FsReview = S.FsModel;
else
    FsReview = Fs;
end

if abs(Fs - FsReview) > 1e-6
    warning('Loaded Fs %.3f differs from Step 6 FsModel %.3f. Using loaded Fs.', Fs, FsReview);
    FsReview = Fs;
end

% Load normalization if possible; otherwise display robust normalized data
[mu, sig] = get_display_mu_sig(S, opts6, size(dataMat,1));

%% ---- Candidate sample mapping ----

hdr = ft_read_header(rawFile);

if isfield(hdr, 'orig') && isfield(hdr.orig, 'sfreq')
    fs_hdr = hdr.orig.sfreq;
else
    fs_hdr = hdr.Fs;
end

try
    first_samp = double(hdr.orig.raw.first_samp);
catch
    first_samp = 0;
end

eventSampleLocal = get_candidate_local_samples(Tfull, FsReview, fs_hdr, first_samp);

valid = eventSampleLocal > 1 & eventSampleLocal < size(dataMat,2);
if any(~valid)
    fprintf('Dropping %d events outside recording boundaries.\n', sum(~valid));
    Tfull = Tfull(valid,:);
    eventSampleLocal = eventSampleLocal(valid);
end

if isempty(eventSampleLocal)
    warning('No valid events to review after boundary check.');
    return;
end

%% ---- Prepare sensor layout ----

dataLite = [];
dataLite.label = chanLabels(megModelIdx);
dataLite.fsample = FsReview;
dataLite.trial = {dataMat(megModelIdx,1:min(size(dataMat,2), round(FsReview)))};
dataLite.time = {(0:size(dataLite.trial{1},2)-1)/FsReview};

layout = get_sensor_layout(dataLite, chanLabels(megModelIdx), opts);

%% ---- Start at first unlabeled event ----

i0 = find(isnan(Tfull.reviewLabel), 1, 'first');

if isempty(i0)
    i0 = 1;
    fprintf('All events already have labels. Starting at event 1.\n');
else
    fprintf('Starting at first unlabeled event: %d/%d\n', i0, height(Tfull));
end

state.i = i0;

%% ---- Create interactive figure ----
% %% ---- Create single-window GUI ------------------------------------------

fig = figure('Name', sprintf('Spike feedback review: %s', step6Base), ...
    'Color','w', ...
    'Units','normalized', ...
    'Position', opts.mainFigPosition, ...
    'WindowKeyPressFcn',@onKey, ...
    'CloseRequestFcn',@onMainClose);

% Panels
controlPanel = uipanel('Parent', fig, ...
    'Units','normalized', ...
    'Position', opts.controlPanelPos, ...
    'BackgroundColor',[0.94 0.94 0.94], ...
    'BorderType','etchedin');

probPanel = uipanel('Parent', fig, ...
    'Units','normalized', ...
    'Position', opts.probPanelPos, ...
    'BackgroundColor','w', ...
    'BorderType','etchedin');

plotPanel = uipanel('Parent', fig, ...
    'Units','normalized', ...
    'Position', opts.plotPanelPos, ...
    'BackgroundColor','w', ...
    'BorderType','none');

navigatorPanel = uipanel('Parent', fig, ...
    'Units','normalized', ...
    'Position', opts.navigatorPanelPos, ...
    'Title','Candidate navigator', ...
    'FontWeight','bold', ...
    'BackgroundColor','w');

% Full-width probability strip
axProb = axes('Parent', probPanel, ...
    'Units','normalized', ...
    'Position', [0.035 0.20 0.94 0.68]);

% Axes inside left plot panel
axTime = axes('Parent', plotPanel, ...
    'Units','normalized', ...
    'Position', opts.axTimePos);

axGfp = axes('Parent', plotPanel, ...
    'Units','normalized', ...
    'Position', opts.axGfpPos);

axMap = axes('Parent', plotPanel, ...
    'Units','normalized', ...
    'Position', opts.axMapPos);

navTable = [];
probAx = axProb;

% The full-width probability strip lives in probPanel. Do not recreate
% axProb inside plotPanel, or the top strip is orphaned/blank.

add_buttons();
createReviewNavigator();

plotCurrent();
updateReviewNavigator();

figure(fig);
% drawnow;
updateReviewNavigator();
drawnow;

uiwait(fig);

%% ===================== Nested callbacks ==============================


    function add_buttons()

        delete(findall(controlPanel, 'Type', 'uicontrol'));

        bg = [0.94 0.94 0.94];

        % ---------------- Row 1: main review controls ----------------
        y1 = 0.58;
        h1 = 0.32;

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Spike (1)', ...
            'Units','normalized','Position',[0.005 y1 0.075 h1], ...
            'Callback',@(src,event)setLabelAndNext(1));

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','No spike (0)', ...
            'Units','normalized','Position',[0.085 y1 0.090 h1], ...
            'Callback',@(src,event)setLabelAndNext(0));

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Unsure (u)', ...
            'Units','normalized','Position',[0.180 y1 0.075 h1], ...
            'Callback',@(src,event)setLabelAndNext(NaN));

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Back (b)', ...
            'Units','normalized','Position',[0.260 y1 0.065 h1], ...
            'Callback',@(src,event)goBack());

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Next (n)', ...
            'Units','normalized','Position',[0.330 y1 0.065 h1], ...
            'Callback',@(src,event)goNext());

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Save (s)', ...
            'Units','normalized','Position',[0.400 y1 0.065 h1], ...
            'Callback',@(src,event)saveReview());

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Save + Quit (q)', ...
            'Units','normalized','Position',[0.470 y1 0.105 h1], ...
            'Callback',@(src,event)saveAndQuit());

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','mbrowse (m)', ...
            'Units','normalized','Position',[0.585 y1 0.105 h1], ...
            'Callback',@(src,event)openMbrowse());

        % ---------------- Row 2: navigation controls ----------------
        y2 = 0.13;
        h2 = 0.30;

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Next unrev (x)', ...
            'Units','normalized','Position',[0.005 y2 0.095 h2], ...
            'Callback',@(src,event)jumpToNextUnreviewed());

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Jump # (j)', ...
            'Units','normalized','Position',[0.105 y2 0.080 h2], ...
            'Callback',@(src,event)jumpToEventNumber());

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Prev spike (z)', ...
            'Units','normalized','Position',[0.190 y2 0.090 h2], ...
            'Callback',@(src,event)jumpToReviewedLabel(1, -1));

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Next spike (c)', ...
            'Units','normalized','Position',[0.285 y2 0.090 h2], ...
            'Callback',@(src,event)jumpToReviewedLabel(1, +1));

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Prev no (a)', ...
            'Units','normalized','Position',[0.380 y2 0.080 h2], ...
            'Callback',@(src,event)jumpToReviewedLabel(0, -1));

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','Next no (d)', ...
            'Units','normalized','Position',[0.465 y2 0.080 h2], ...
            'Callback',@(src,event)jumpToReviewedLabel(0, +1));

        uicontrol('Parent',controlPanel,'Style','pushbutton','String','High p unrev (h)', ...
            'Units','normalized','Position',[0.555 y2 0.10 h2], ...
            'Callback',@(src,event)jumpToHighestProbUnreviewed());

        % ---------------- Help text, right side ----------------
        uicontrol('Parent',controlPanel,'Style','text', ...
            'String','Keys: 1 spike | 0 no-spike | u unsure | n next | b back | x unreviewed | h high-p | j jump | m mbrowse | s save | q quit', ...
            'Units','normalized','Position',[0.665 y2 0.330 h2], ...
            'BackgroundColor',bg, ...
            'HorizontalAlignment','left', ...
            'FontSize',8);
    end

    function plotCurrent()

        if state.i < 1
            state.i = 1;
        end

        if state.i > height(Tfull)
            state.i = height(Tfull);
        end

        c = eventSampleLocal(state.i);

        nPre  = round(abs(opts.reviewWinSec(1)) * FsReview);
        nPost = round(abs(opts.reviewWinSec(2)) * FsReview);

        idx = (c-nPre):(c+nPost);

        if min(idx) < 1 || max(idx) > size(dataMat,2)
            warning('Candidate window outside range.');
            return;
        end

        t = ((idx(:) - c) ./ FsReview)';

        %% --------------------------------------------------------------------
        % MEG-only data for GFP and sensor-layout waveform map
        % ---------------------------------------------------------------------

        % Safety fallback if megModelIdx was not defined outside
        if ~exist('megModelIdx', 'var') || isempty(megModelIdx)
            chInfoLocal = get_channel_groups(chanLabels);
            megModelIdxLocal = chInfoLocal.megIdx;
            megModelIdxLocal = megModelIdxLocal(1:min(306, numel(megModelIdxLocal)));
        else
            megModelIdxLocal = megModelIdx;
        end

        XmegRaw = single(dataMat(megModelIdxLocal, idx));  % MEG only, Cmeg x T

        % Normalize MEG channels using model mu/sig if dimensions match
        if numel(mu) == size(XmegRaw,1)
            XmegZ = bsxfun(@rdivide, bsxfun(@minus, XmegRaw, mu), sig);
        else
            XmegZ = robust_display_scale(XmegRaw, 'per_channel', 98, [], 1.0);
        end

        %% --------------------------------------------------------------------
        % Time-course panel: selected sensor type, e.g. MAG + ECG
        % ---------------------------------------------------------------------

        tcIdx = select_timecourse_channels(chanLabels, opts.timeCourseSensorType);

        if isempty(tcIdx)
            warning('No channels found for timeCourseSensorType=%s. Falling back to MEG.', ...
                opts.timeCourseSensorType);
            tcIdx = megModelIdxLocal;
        end

        % Limit number of traces shown
        if isfield(opts, 'maxTimeCourseChannels') && numel(tcIdx) > opts.maxTimeCourseChannels
            tcIdx = tcIdx(1:opts.maxTimeCourseChannels);
        end

        XtcRaw = single(dataMat(tcIdx, idx));
        tcLabels = chanLabels(tcIdx);

        XtcDisp = robust_display_scale( ...
            XtcRaw, ...
            opts.timeCourseScaleMode, ...
            opts.timeCourseScalePrctile, ...
            opts.timeCourseManualScale, ...
            opts.timeCourseAmpMultiplier);

        %% --------------------------------------------------------------------
        % GFP from MEG only
        % ---------------------------------------------------------------------

        gfp = sqrt(mean(XmegZ.^2, 1));

        % Local GFP peak around event, mainly for visual reference
        peakLo = max(1, round((opts.peakSearchSec(1) - opts.reviewWinSec(1)) * FsReview));
        peakHi = min(numel(t), round((opts.peakSearchSec(2) - opts.reviewWinSec(1)) * FsReview));

        if peakHi > peakLo
            [~, iLocal] = max(gfp(peakLo:peakHi));
            iPeak = peakLo + iLocal - 1;
        else
            iPeak = nearest_idx(t, 0);
        end

        iZero = nearest_idx(t, 0);

        %% --------------------------------------------------------------------
        % Figure layout
        cla(axTime);
        % cla(axProb);
        cla(axGfp);
        cla(axMap);

        probTxt = '';
        if ismember('probSpike', Tfull.Properties.VariableNames)
            probTxt = sprintf(' | p=%.3f', Tfull.probSpike(state.i));
        end

        labelTxt = label_to_text(Tfull.reviewLabel(state.i));

        %% --------------------------------------------------------------------
        % Large selected time-course panel
        % ---------------------------------------------------------------------

        plot_selected_timecourses( ...
            axTime, t, XtcDisp, tcLabels, ...
            state.i, height(Tfull), probTxt, labelTxt, opts);
        updateProbabilityNavigator();

        %% --------------------------------------------------------------------
        % GFP panel
        % ---------------------------------------------------------------------

        axes(axGfp);
        cla(axGfp);

        plot(t, gfp, 'k', 'LineWidth', 1.2);
        hold on;
        xline(0, 'g', 'LineWidth', 1.5);              % candidate/spike time
        xline(t(iPeak), 'r--', 'LineWidth', 1.2);     % local GFP peak
        scatter(t(iPeak), gfp(iPeak), 50, 'r', 'filled');
        scatter(t(iZero), gfp(iZero), 35, 'g', 'filled');
        hold off;

        xlabel('Time around event (s)');
        ylabel('GFP / normalized energy');
        % title('Global field power');
        title('GFP');
        set(axGfp, 'FontSize', 7);
        grid on;

        %% --------------------------------------------------------------------
        % Sensor waveform map: only -0.15 to +0.15 sec
        % ---------------------------------------------------------------------
        % Sensor waveform map: only selected MEG sensor type and only -0.15 to +0.15 sec
        mapMask = t >= opts.sensorMapWinSec(1) & t <= opts.sensorMapWinSec(2);

        if ~any(mapMask)
            warning('sensorMapWinSec has no samples. Using full time window for sensor map.');
            tMapFull = t;
            XmapFull = XmegZ;
        else
            tMapFull = t(mapMask);
            XmapFull = XmegZ(:, mapMask);
        end

        % Select which MEG channels to show in the sensor map
        mapIdx = select_sensor_map_channels(layout.label, opts.sensorMapSensorType);

        if isempty(mapIdx)
            warning('No channels found for sensorMapSensorType=%s. Using all MEG channels.', ...
                opts.sensorMapSensorType);
            mapIdx = 1:size(XmapFull,1);
        end

        Xmap = XmapFull(mapIdx,:);
        layoutMap = layout;
        layoutMap.label = layout.label(mapIdx);
        layoutMap.pos   = layout.pos(mapIdx,:);

        plot_sensor_waveform_map_mbrowse(axMap, tMapFull, Xmap, layoutMap, opts);

        title(axMap, sprintf('Sensor waveform map %.2f to %.2f s', ...
            opts.sensorMapWinSec(1), opts.sensorMapWinSec(2)));

        % drawnow;
        updateProbabilityNavigator();
        updateReviewNavigator();
        drawnow;

    end
    function createReviewNavigator()

        delete(findall(navigatorPanel, 'Type', 'uitable'));
        delete(findall(navigatorPanel, 'Type', 'axes'));
        delete(findall(navigatorPanel, 'Type', 'uicontrol'));

        % Probability bar plot, top of navigator panel
        % if isfield(opts, 'showProbabilityNavigator') && opts.showProbabilityNavigator
        %     probAx = axes('Parent', navigatorPanel, ...
        %         'Units','normalized', ...
        %         'Position', opts.probNavigatorAxesPosition);
        % end

        % Instruction text
        uicontrol('Parent', navigatorPanel, ...
            'Style','text', ...
            'String','Click table row or probability bar to jump', ...
            'Units','normalized', ...
            'Position',[0.03 0.625 0.94 0.035], ...
            'BackgroundColor','w', ...
            'FontWeight','bold', ...
            'HorizontalAlignment','left');

        % Candidate table
        % navTable = uitable('Parent', navigatorPanel, ...
        %     'Units','normalized', ...
        %     'Position', opts.navigatorTablePosition, ...
        %     'CellSelectionCallback', @onNavigatorSelect);

        % navTable = uitable('Parent', navigatorPanel, ...
        %     'Units','normalized', ...
        %     'Position', opts.navigatorTablePosition, ...
        %     'FontSize', opts.navigatorTableFontSize, ...
        %     'CellSelectionCallback', @onNavigatorSelect);

        navTable = uitable('Parent', navigatorPanel, ...
            'Units','normalized', ...
            'Position', opts.navigatorTablePosition, ...
            'FontSize', opts.navigatorTableFontSize, ...
            'RowName', [], ...
            'CellSelectionCallback', @onNavigatorSelect);
    end

    function updateProbabilityNavigator()

        if isempty(probAx) || ~isvalid(probAx)
            return;
        end

        if ~ismember('probSpike', Tfull.Properties.VariableNames)
            cla(probAx);
            title(probAx, 'No probSpike');
            return;
        end

        N = height(Tfull);
        x = 1:N;
        p = double(Tfull.probSpike(:));

        cla(probAx, 'reset');

        hBar = bar(probAx, x, p, opts.probNavigatorBarWidth, ...
            'FaceColor','flat', ...
            'EdgeColor', opts.probNavigatorEdgeColor, ...
            'LineWidth', opts.probNavigatorLineWidth);

        % Color bars by review status:
        % gray = unreviewed, red = reviewed spike, blue = reviewed no-spike,
        % green = current candidate
        C = repmat([0.70 0.70 0.70], N, 1);

        idxSpike = Tfull.reviewLabel == 1;
        idxNo    = Tfull.reviewLabel == 0;
        idxUnrev = isnan(Tfull.reviewLabel);

        C(idxUnrev,:) = repmat([0.70 0.70 0.70], sum(idxUnrev), 1);
        C(idxSpike,:) = repmat([0.85 0.15 0.15], sum(idxSpike), 1);
        C(idxNo,:)    = repmat([0.15 0.35 0.85], sum(idxNo), 1);

        if state.i >= 1 && state.i <= N
            C(state.i,:) = [0.00 0.65 0.00];
        end

        hBar.CData = C;

        set(hBar, 'PickableParts','all', 'HitTest','on', ...
            'ButtonDownFcn', @onProbabilityBarClick);

        set(probAx, 'ButtonDownFcn', @onProbabilityBarClick);

        hold(probAx, 'on');

        xline(probAx, state.i, 'g-', 'LineWidth', 1.4);

        if isfield(opts, 'probNavigatorShowThreshold') && opts.probNavigatorShowThreshold
            if isfield(opts6, 'probThreshold')
                yline(probAx, opts6.probThreshold, 'k--', ...
                    sprintf('thr %.2f', opts6.probThreshold), ...
                    'FontSize', 6, 'LineWidth', 0.8);
            end
        end

        hold(probAx, 'off');

        xlim(probAx, [0.5 N+0.5]);
        ylim(probAx, [0 1]);

        if ~isfield(opts, 'probNavigatorXLabelEvery')
            opts.probNavigatorXLabelEvery = 1;
        end

        xTicks = 1:opts.probNavigatorXLabelEvery:N;

        xticks(probAx, xTicks);
        xticklabels(probAx, string(xTicks));

        if isfield(opts, 'probNavigatorXTickRotation')
            probAx.XTickLabelRotation = opts.probNavigatorXTickRotation;
        end

        xlabel(probAx, 'Candidate #', 'FontSize', 6);
        ylabel(probAx, 'p', 'FontSize', 6);

        title(probAx, sprintf('Spike probability | current=%d/%d', state.i, N), ...
            'FontSize', 7, 'FontWeight','bold');

        set(probAx, ...
            'FontSize', 6, ...
            'YTick', [0 0.5 1], ...
            'Box', 'on', ...
            'TickDir','out');

        grid(probAx, 'on');
    end

    function updateReviewNavigator()

        if isempty(navTable) || ~isvalid(navTable)
            return;
        end

        N = height(Tfull);

        idxCol = (1:N)';

        timeCol = nan(N,1);
        if ismember('eventTimeSec', Tfull.Properties.VariableNames)
            timeCol = Tfull.eventTimeSec;
        elseif ismember('modelCenterTimeSec', Tfull.Properties.VariableNames)
            timeCol = Tfull.modelCenterTimeSec;
        end

        probCol = nan(N,1);
        if ismember('probSpike', Tfull.Properties.VariableNames)
            probCol = Tfull.probSpike;
        end

        labelCol = strings(N,1);
        for ii = 1:N
            if isnan(Tfull.reviewLabel(ii))
                labelCol(ii) = "unreviewed";
            elseif Tfull.reviewLabel(ii) == 1
                labelCol(ii) = "SPIKE";
            elseif Tfull.reviewLabel(ii) == 0
                labelCol(ii) = "NO-SPIKE";
            else
                labelCol(ii) = "unknown";
            end
        end

        currentCol = repmat("", N, 1);
        currentCol(state.i) = "<-- current";

        D = cell(N,5);
        for ii = 1:N
            D{ii,1} = idxCol(ii);
            D{ii,2} = timeCol(ii);
            D{ii,3} = probCol(ii);
            D{ii,4} = char(labelCol(ii));
            D{ii,5} = char(currentCol(ii));
        end

        % set(navTable, ...
        %     'Data', D, ...
        %     'ColumnName', {'Idx','Time_sec','ProbSpike','Review','Current'}, ...
        %     'ColumnWidth', {35,65,65,75,75});

        set(navTable, ...
            'Data', D, ...
            'ColumnName', {'Idx','Time_sec','ProbSpike','Review','Current'}, ...
            'ColumnWidth', opts.navigatorColumnWidth, ...
            'FontSize', opts.navigatorTableFontSize);

        navTable.UserData.rowIdx = idxCol;

        if isfield(opts, 'navigatorTableFontSize')
            set(navTable, 'FontSize', opts.navigatorTableFontSize);
        else
            set(navTable, 'FontSize', 7);
        end

        navTable.UserData.rowIdx = idxCol;
        updateProbabilityNavigator();
    end

    function onNavigatorSelect(src, event)

        if isempty(event.Indices)
            return;
        end

        row = event.Indices(1);

        if isfield(src.UserData, 'rowIdx')
            newIdx = src.UserData.rowIdx(row);
        else
            newIdx = row;
        end

        if newIdx >= 1 && newIdx <= height(Tfull)
            state.i = newIdx;
            figure(fig);
            plotCurrent();
            updateReviewNavigator();
        end
    end

    function jumpToEventNumber()

        answer = inputdlg( ...
            sprintf('Enter event number, 1 to %d:', height(Tfull)), ...
            'Jump to event', ...
            [1 40], ...
            {num2str(state.i)});

        if isempty(answer)
            return;
        end

        newIdx = str2double(answer{1});

        if isnan(newIdx) || newIdx < 1 || newIdx > height(Tfull)
            warning('Invalid event number.');
            return;
        end

        state.i = round(newIdx);
        plotCurrent();
        updateReviewNavigator();
    end

    function jumpToNextUnreviewed()

        N = height(Tfull);
        allIdx = (1:N)';

        idx = find(isnan(Tfull.reviewLabel) & allIdx > state.i, 1, 'first');

        % wrap around
        if isempty(idx)
            idx = find(isnan(Tfull.reviewLabel), 1, 'first');
        end

        if isempty(idx)
            fprintf('No unreviewed events remain.\n');
            return;
        end

        state.i = idx;
        plotCurrent();
        updateReviewNavigator();
    end

    function jumpToReviewedLabel(labelValue, direction)

        N = height(Tfull);
        allIdx = (1:N)';

        if direction > 0
            idx = find(Tfull.reviewLabel == labelValue & allIdx > state.i, 1, 'first');

            % wrap around
            if isempty(idx)
                idx = find(Tfull.reviewLabel == labelValue, 1, 'first');
            end
        else
            idx = find(Tfull.reviewLabel == labelValue & allIdx < state.i, 1, 'last');

            % wrap around
            if isempty(idx)
                idx = find(Tfull.reviewLabel == labelValue, 1, 'last');
            end
        end

        if isempty(idx)
            if labelValue == 1
                fprintf('No reviewed spike events yet.\n');
            else
                fprintf('No reviewed no-spike events yet.\n');
            end
            return;
        end

        state.i = idx;
        plotCurrent();
        updateReviewNavigator();
    end

    function setLabelAndNext(labelValue)

        labeledIdx = state.i;

        if isnan(labelValue)
            Tfull.reviewLabel(labeledIdx) = NaN;
            Tfull.reviewType{labeledIdx} = 'unsure';
            Tfull.reviewNotes{labeledIdx} = 'unsure';
        elseif labelValue == 1
            Tfull.reviewLabel(labeledIdx) = 1;
            Tfull.reviewType{labeledIdx} = 'true_spike';
            Tfull.reviewNotes{labeledIdx} = '';
        elseif labelValue == 0
            Tfull.reviewLabel(labeledIdx) = 0;
            Tfull.reviewType{labeledIdx} = 'false_positive';
            Tfull.reviewNotes{labeledIdx} = '';
        end

        fprintf('Event %d labeled: %s\n', labeledIdx, label_to_text(labelValue));

        if state.i < height(Tfull)
            state.i = state.i + 1;
        end

        plotCurrent();
        updateReviewNavigator();
        drawnow;

    end

    function goBack()
        state.i = max(1, state.i - 1);
        plotCurrent();
        updateReviewNavigator();
        drawnow;
    end

    function goNext()
        state.i = min(height(Tfull), state.i + 1);
        plotCurrent();
        updateReviewNavigator();
        drawnow;
    end

    function saveReview()

        Treviewed = Tfull;

        save(outMat, ...
            'Treviewed','eventSampleLocal','opts','opts6','rawFile','FsReview','chanLabels', ...
            '-v7.3');

        writetable(Treviewed, outCsv);
        writetable(Treviewed, outXlsx);

        if isfield(opts, 'exportReviewedEve') && opts.exportReviewedEve

            reviewedEveFile = export_reviewed_eve_file( ...
                Treviewed, eventSampleLocal, rawFile, FsReview, ...
                opts.reviewedSpikeCode, opts.reviewedNoSpikeCode, outMat);

            fprintf('Saved reviewed EVE file:\n%s\n', reviewedEveFile);
        end

        nReviewed = sum(~isnan(Treviewed.reviewLabel));
        nSpike = sum(Treviewed.reviewLabel == 1);
        nNo = sum(Treviewed.reviewLabel == 0);

        fprintf('\nSaved review:\n%s\n%s\n%s\n', outMat, outCsv, outXlsx);
        fprintf('Reviewed: %d/%d | Spike=%d | No-spike=%d\n', ...
            nReviewed, height(Treviewed), nSpike, nNo);
    end

    function openMbrowse()

        % Save current review first, including reviewed .mat/.csv/.xlsx
        saveReview();

        fprintf('\nLaunching mbrowse for raw FIF:\n%s\n', rawFile);

        % If you are also exporting reviewed EVE files, print the likely file path
        [outDirTmp, outBaseTmp, ~] = fileparts(outMat);
        reviewedEveFile = fullfile(outDirTmp, [outBaseTmp '_reviewed_1111_2222.eve']);

        if isfile(reviewedEveFile)
            fprintf('Reviewed EVE file available:\n%s\n', reviewedEveFile);
        else
            fprintf('Reviewed EVE file not found yet. If needed, save/export reviewed EVE first.\n');
        end

        % Launch mbrowse in background so MATLAB stays interactive.
        % Correct command is "mbrowse", not "mbrowose".
        cmd = sprintf('mbrowse "%s" &', rawFile);

        [status, msg] = system(cmd);

        if status ~= 0
            warning('Could not launch mbrowse.\nCommand: %s\nMessage:\n%s', cmd, msg);
        else
            fprintf('mbrowse launched.\n');
        end
    end

    function saveAndQuit()

        saveReview();

        if exist('fig','var') && ~isempty(fig) && isvalid(fig)
            uiresume(fig);
            delete(fig);
        end
    end


    function onKey(~, event)

        switch lower(event.Key)

            case {'1','numpad1'}
                setLabelAndNext(1);

            case {'0','numpad0'}
                setLabelAndNext(0);

            case {'u'}
                setLabelAndNext(NaN);

            case {'n','rightarrow','space'}
                goNext();

            case {'b','leftarrow'}
                goBack();

            case {'s'}
                saveReview();

            case {'q','escape'}
                saveAndQuit();

            case {'j'}
                jumpToEventNumber();

            case {'x'}
                jumpToNextUnreviewed();

            case {'z'}
                jumpToReviewedLabel(1, -1);   % previous spike

            case {'c'}
                jumpToReviewedLabel(1, +1);   % next spike

            case {'a'}
                jumpToReviewedLabel(0, -1);   % previous no-spike

            case {'d'}
                jumpToReviewedLabel(0, +1);   % next no-spike
            case {'m'}
                openMbrowse();
            case {'h'}
                jumpToHighestProbUnreviewed();
        end
    end
end

%% ========================================================================
% Helper functions
% ========================================================================

function [dataMat, chanLabels, Fs] = load_meg_for_review(rawFile, opts, opts6)

cfg = [];
cfg.dataset = rawFile;
cfg.continuous = 'yes';

if isfield(opts, 'reviewLoadChannels')
    cfg.channel = opts.reviewLoadChannels;
else
    cfg.channel = 'MEG';
end

% Match Step 6 display preprocessing if available
if opts.matchStep6Bandpass && isfield(opts6, 'useBandpass') && opts6.useBandpass

    fprintf('Applying display bandpass %.1f-%.1f Hz...\n', ...
        opts6.bpFreq(1), opts6.bpFreq(2));

    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = opts6.bpFreq;
    cfg.bpfilttype = 'but';

    if isfield(opts6, 'bpOrder')
        cfg.bpfiltord = opts6.bpOrder;
    else
        cfg.bpfiltord = 4;
    end

    cfg.bpfiltdir = 'twopass';
end

try
    data = ft_preprocessing(cfg);
catch ME
    warning('Could not load requested review channels. Falling back to MEG only.\n%s', ME.message);

    cfg = [];
    cfg.dataset = rawFile;
    cfg.continuous = 'yes';
    cfg.channel = 'MEG';

    if opts.matchStep6Bandpass && isfield(opts6, 'useBandpass') && opts6.useBandpass
        cfg.bpfilter   = 'yes';
        cfg.bpfreq     = opts6.bpFreq;
        cfg.bpfilttype = 'but';
        cfg.bpfiltord  = opts6.bpOrder;
        cfg.bpfiltdir  = 'twopass';
    end

    data = ft_preprocessing(cfg);
end

Fs = data.fsample;
dataMat = data.trial{1};
chanLabels = data.label(:);
end

function eventSampleLocal = get_candidate_local_samples(Tfull, Fs, fs_hdr, first_samp)

names = Tfull.Properties.VariableNames;

if ismember('refinedSampleLocal', names)
    eventSampleLocal = round(double(Tfull.refinedSampleLocal(:)));

elseif ismember('modelCenterSampleLocal', names)
    eventSampleLocal = round(double(Tfull.modelCenterSampleLocal(:)));

elseif ismember('eventSample', names)
    eventSampleAbs = double(Tfull.eventSample(:));
    eventSampleLocal = round((eventSampleAbs - first_samp) ./ fs_hdr .* Fs) + 1;

elseif ismember('sample', names)
    eventSampleAbs = double(Tfull.sample(:));
    eventSampleLocal = round((eventSampleAbs - first_samp) ./ fs_hdr .* Fs) + 1;

elseif ismember('eventTimeSec', names)
    eventTimeAbs = double(Tfull.eventTimeSec(:));
    eventTimeRel = eventTimeAbs - first_samp ./ fs_hdr;
    eventSampleLocal = round(eventTimeRel .* Fs) + 1;

elseif ismember('time_sec', names)
    eventTimeAbs = double(Tfull.time_sec(:));
    eventTimeRel = eventTimeAbs - first_samp ./ fs_hdr;
    eventSampleLocal = round(eventTimeRel .* Fs) + 1;

else
    error('Could not find event sample/time columns in Tfull.');
end
end

function [mu, sig] = get_display_mu_sig(S, opts6, C)

mu = [];
sig = [];

if isfield(S, 'opts') && isfield(opts6, 'datasetFile') && isfile(opts6.datasetFile)
    [mu, sig] = load_mu_sig(opts6.datasetFile);
end

if isempty(mu) || isempty(sig) || numel(mu) ~= C
    mu = zeros(C,1,'single');
    sig = ones(C,1,'single');
else
    mu = single(mu(:));
    sig = single(sig(:));
end

sig(sig == 0 | isnan(sig)) = 1;
mu(isnan(mu)) = 0;
end

function [mu, sig] = load_mu_sig(datasetFile)

info = whos('-file', datasetFile);
names = {info.name};

if all(ismember({'mu','sig'}, names))
    S = load(datasetFile, 'mu','sig');
    mu = S.mu;
    sig = S.sig;

elseif ismember('dataset', names)
    S = load(datasetFile, 'dataset');
    if isfield(S.dataset,'mu') && isfield(S.dataset,'sig')
        mu = S.dataset.mu;
        sig = S.dataset.sig;
    else
        mu = [];
        sig = [];
    end
else
    mu = [];
    sig = [];
end
end

function txt = label_to_text(v)

if isnan(v)
    txt = 'UNSURE';
elseif v == 1
    txt = 'SPIKE';
elseif v == 0
    txt = 'NO-SPIKE';
else
    txt = 'UNKNOWN';
end
end

function layout = get_sensor_layout(data, chanLabels, opts)

layout = struct();

try
    cfg = [];
    cfg.layout = opts.layoutFile;
    cfg.channel = chanLabels;
    cfg.skipscale = 'yes';
    cfg.skipcomnt = 'yes';

    lay = ft_prepare_layout(cfg, data);

    labels = lay.label(:);
    pos = lay.pos;

    bad = strcmp(labels,'COMNT') | strcmp(labels,'SCALE');
    labels = labels(~bad);
    pos = pos(~bad,:);

    [tf, loc] = ismember(chanLabels, labels);

    layout.label = chanLabels;
    layout.pos = nan(numel(chanLabels),2);
    layout.pos(tf,:) = pos(loc(tf),:);

    miss = ~tf;
    if any(miss)
        warning('%d channels missing from layout. Filling missing channels on circle.', sum(miss));

        theta = linspace(0,2*pi,sum(miss)+1)';
        theta(end) = [];
        layout.pos(miss,:) = [cos(theta), sin(theta)];
    end

catch ME
    warning('Could not prepare FieldTrip layout %s: %s. Using circle fallback.', ...
        opts.layoutFile, ME.message);

    n = numel(chanLabels);
    theta = linspace(0,2*pi,n+1)';
    theta(end) = [];

    layout.label = chanLabels;
    layout.pos = [cos(theta), sin(theta)];
end

p = layout.pos;
p = p - mean(p,1,'omitnan');

mx = max(abs(p), [], 'all');

if mx > 0
    p = p ./ mx;
end

layout.pos = p;
end

function plot_sensor_waveform_map_mbrowse(ax, t, Xz, layout, opts)

if ~isfield(opts, 'sensorMapTimeScale')
    opts.sensorMapTimeScale = 0.12;
end

if ~isfield(opts, 'sensorMapAmpScale')
    opts.sensorMapAmpScale = 0.045;
end

if ~isfield(opts, 'maxSensorWaveforms')
    opts.maxSensorWaveforms = size(Xz,1);
end

if ~isfield(opts, 'sensorMapShowLabels')
    opts.sensorMapShowLabels = true;
end

if ~isfield(opts, 'sensorMapLabelEvery')
    opts.sensorMapLabelEvery = 6;
end

if ~isfield(opts, 'sensorMapLabelFontSize')
    opts.sensorMapLabelFontSize = 5;
end

axes(ax);
cla(ax);

C = size(Xz,1);
% maxSensors = min(opts.maxSensorWaveforms, C);

maxSensors = min(opts.maxSensorWaveforms, C);
if C <= opts.maxSensorWaveforms
    chUse = (1:C)';
else
    score = max(abs(Xz), [], 2);
    [~, ord] = sort(score, 'descend');
    chUse = sort(ord(1:maxSensors));
end

% score = max(abs(Xz), [], 2);
% [~, ord] = sort(score, 'descend');
% chUse = sort(ord(1:maxSensors));

pos = layout.pos(chUse,:);
Xuse = Xz(chUse,:);

ampScale = prctile(abs(Xuse(:)), 98);
if ampScale == 0 || isnan(ampScale)
    ampScale = 1;
end

tNorm = (t - min(t)) ./ (max(t) - min(t));
tNorm = (tNorm - 0.5) * opts.sensorMapTimeScale;

iZero = nearest_idx(t, 0);

hold on;

for ii = 1:numel(chUse)

    x0 = pos(ii,1);
    y0 = pos(ii,2);

    y = Xuse(ii,:);
    y = y - median(y, 'omitnan');
    y = y ./ ampScale;
    y = y * opts.sensorMapAmpScale;

    % waveform
    plot(x0 + tNorm, y0 + y, 'b', 'LineWidth', 0.5);

    % spike-time marker at t = 0
    xCenter = x0 + tNorm(iZero);
    plot([xCenter xCenter], ...
        [y0 - opts.sensorMapAmpScale, y0 + opts.sensorMapAmpScale], ...
        'g', 'LineWidth', 0.5);

    % optional sensor labels
    if opts.sensorMapShowLabels && mod(ii-1, opts.sensorMapLabelEvery) == 0
        if isfield(layout, 'label') && numel(layout.label) >= chUse(ii)
            lab = layout.label{chUse(ii)};
        else
            lab = sprintf('%d', chUse(ii));
        end

        text(x0, y0 + 1.6*opts.sensorMapAmpScale, lab, ...
            'FontSize', opts.sensorMapLabelFontSize, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', ...
            'Color','k', ...
            'Interpreter','none');
    end
end

hold off;

axis equal off;
end

function chInfo = get_channel_groups(chanLabels)

labels = upper(string(chanLabels));

isMEG = startsWith(labels, "MEG");

% Neuromag convention:
% MEGxxxx1 = magnetometer
% MEGxxxx2/3 = planar gradiometers
isMag  = isMEG & endsWith(labels, "1");
isGrad = isMEG & (endsWith(labels, "2") | endsWith(labels, "3"));

isECG = contains(labels, "ECG");
isEOG = contains(labels, "EOG");

chInfo.megIdx  = find(isMEG);
chInfo.magIdx  = find(isMag);
chInfo.gradIdx = find(isGrad);
chInfo.ecgIdx  = find(isECG);
chInfo.eogIdx  = find(isEOG);
end

function idx = select_timecourse_channels(chanLabels, sensorType)

ch = get_channel_groups(chanLabels);

switch lower(strrep(sensorType, '+', '_'))

    case {'mag','magnetometer','magnetometers'}
        idx = ch.magIdx;

    case {'grad','gradiometer','gradiometers'}
        idx = ch.gradIdx;

    case {'meg','allmeg','all_meg'}
        idx = ch.megIdx;

    case {'mag_ecg','magandecg'}
        idx = [ch.magIdx; ch.ecgIdx];

    case {'meg_ecg','megandecg'}
        idx = [ch.megIdx; ch.ecgIdx];

    case {'ecg'}
        idx = ch.ecgIdx;

    case {'eog'}
        idx = ch.eogIdx;

    case {'all'}
        idx = (1:numel(chanLabels))';

    otherwise
        warning('Unknown sensor type: %s. Using MAG.', sensorType);
        idx = ch.magIdx;
end

idx = unique(idx(:), 'stable');
end

function Xdisp = robust_display_scale(Xraw, scaleMode, pct, manualScale, ampMultiplier)

Xraw = single(Xraw);

if nargin < 2 || isempty(scaleMode)
    scaleMode = 'per_channel';
end

if nargin < 3 || isempty(pct)
    pct = 98;
end

if nargin < 4
    manualScale = [];
end

if nargin < 5 || isempty(ampMultiplier)
    ampMultiplier = 1.0;
end

X = Xraw;

% Remove median per channel for display
X = bsxfun(@minus, X, median(X, 2, 'omitnan'));

switch lower(scaleMode)

    case 'per_channel'

        scaleVal = prctile(abs(X), pct, 2);
        scaleVal(scaleVal == 0 | isnan(scaleVal)) = 1;
        Xdisp = bsxfun(@rdivide, X, scaleVal);

    case 'global'

        scaleVal = prctile(abs(X(:)), pct);
        if scaleVal == 0 || isnan(scaleVal)
            scaleVal = 1;
        end
        Xdisp = X ./ scaleVal;

    case 'manual'

        if isempty(manualScale) || manualScale == 0
            error('manualScale must be non-empty when scaleMode = manual.');
        end
        Xdisp = X ./ manualScale;

    otherwise

        warning('Unknown scaleMode=%s. Using per_channel.', scaleMode);
        scaleVal = prctile(abs(X), pct, 2);
        scaleVal(scaleVal == 0 | isnan(scaleVal)) = 1;
        Xdisp = bsxfun(@rdivide, X, scaleVal);
end

Xdisp = Xdisp .* ampMultiplier;
end

function plot_selected_timecourses(ax, t, Xdisp, labels, iEvent, nEvents, probTxt, labelTxt, opts)

axes(ax);
cla(ax);

nShow = size(Xdisp,1);
offset = opts.timeCourseOffset;

yBase = zeros(nShow,1);

hold on;

for kk = 1:nShow

    yBase(kk) = offset*(nShow-kk);
    y = Xdisp(kk,:) + yBase(kk);

    lab = upper(string(labels{kk}));

    if contains(lab, "ECG")
        plot(t, y, 'r', 'LineWidth', 1.0);
    else
        plot(t, y, 'b', 'LineWidth', 0.6);
    end
end

xline(0, 'g', 'LineWidth', 1.5);

hold off;

xlabel('Time around event (s)');
ylabel(sprintf('%s channels', opts.timeCourseSensorType));

title(sprintf('Event %d/%d%s | label=%s | %s', ...
    iEvent, nEvents, probTxt, labelTxt, opts.timeCourseSensorType), ...
    'Interpreter','none');

grid on;
axis tight;



% Add y-axis sensor labels
if isfield(opts, 'timeCourseShowYLabels') && opts.timeCourseShowYLabels

    if ~isfield(opts, 'timeCourseYLabelEvery')
        opts.timeCourseYLabelEvery = 5;
    end

    if ~isfield(opts, 'timeCourseLabelFontSize')
        opts.timeCourseLabelFontSize = 6;
    end

    tickIdx = 1:opts.timeCourseYLabelEvery:nShow;

    % Always label ECG channels if present
    isEcg = contains(upper(string(labels)), "ECG");
    tickIdx = unique([tickIdx(:); find(isEcg(:))], 'stable');

    % MATLAB yticks must be increasing.
    tickVals = yBase(tickIdx);
    tickLabs = labels(tickIdx);

    [tickValsSorted, ordTick] = sort(tickVals, 'ascend');
    tickLabsSorted = tickLabs(ordTick);

    yticks(ax, tickValsSorted);
    yticklabels(ax, tickLabsSorted);

    set(ax, 'FontSize', opts.timeCourseLabelFontSize);
end
end

function idx = nearest_idx(t, val)

[~, idx] = min(abs(t - val));
end


function idx = select_sensor_map_channels(chanLabels, sensorType)

labels = upper(string(chanLabels));
isMEG = startsWith(labels, "MEG");

% Neuromag convention:
% MEGxxxx1 = magnetometer
% MEGxxxx2/3 = planar gradiometers
isMag  = isMEG & endsWith(labels, "1");
isGrad = isMEG & (endsWith(labels, "2") | endsWith(labels, "3"));

switch lower(sensorType)

    case {'mag','magnetometer','magnetometers'}
        idx = find(isMag);

    case {'grad','gradiometer','gradiometers'}
        idx = find(isGrad);

    case {'meg','allmeg','all_meg'}
        idx = find(isMEG);

    otherwise
        warning('Unknown sensorMapSensorType=%s. Using MAG.', sensorType);
        idx = find(isMag);
end

idx = idx(:);
end

function reviewedEveFile = export_reviewed_eve_file( ...
    Treviewed, eventSampleLocal, rawFile, FsReview, spikeCode, noSpikeCode, outMat)

[outDir, outBase, ~] = fileparts(outMat);

reviewedEveFile = fullfile(outDir, [outBase '_reviewed_1111_2222.eve']);

hdr = ft_read_header(rawFile);

if isfield(hdr, 'orig') && isfield(hdr.orig, 'sfreq')
    fs_hdr = hdr.orig.sfreq;
else
    fs_hdr = hdr.Fs;
end

try
    first_samp = double(hdr.orig.raw.first_samp);
catch
    first_samp = 0;
end

fid = fopen(reviewedEveFile, 'w');
assert(fid > 0, 'Cannot open reviewed EVE file: %s', reviewedEveFile);

% Header/seed line for mbrowse compatibility
fprintf(fid, '%d\t%f\t%d\t%d\t%s\n', ...
    first_samp, first_samp/fs_hdr, 0, 0, 'reviewed');

nSpike = 0;
nNoSpike = 0;

for ii = 1:height(Treviewed)

    lab = Treviewed.reviewLabel(ii);

    % Skip unsure/unreviewed rows
    if isnan(lab)
        continue;
    end

    if lab == 1
        code = spikeCode;
        nSpike = nSpike + 1;
    elseif lab == 0
        code = noSpikeCode;
        nNoSpike = nNoSpike + 1;
    else
        continue;
    end

    % Prefer absolute sample/time from Step 6 table if present
    if ismember('eventSample', Treviewed.Properties.VariableNames)
        samp = round(Treviewed.eventSample(ii));
    elseif ismember('sample', Treviewed.Properties.VariableNames)
        samp = round(Treviewed.sample(ii));
    else
        samp = round((eventSampleLocal(ii)-1) ./ FsReview .* fs_hdr) + first_samp;
    end

    if ismember('eventTimeSec', Treviewed.Properties.VariableNames)
        tsec = Treviewed.eventTimeSec(ii);
    elseif ismember('time_sec', Treviewed.Properties.VariableNames)
        tsec = Treviewed.time_sec(ii);
    else
        tRel = double(eventSampleLocal(ii)-1) ./ FsReview;
        tsec = tRel + first_samp ./ fs_hdr;
    end

    fprintf(fid, '%d\t%f\t%d\t%d\n', ...
        round(samp), tsec, 0, code);
end

fclose(fid);

fprintf('Reviewed EVE contents: spike=%d, no-spike=%d\n', nSpike, nNoSpike);
end

function onMainClose(~,~)

saveReview();

if exist('fig','var') && ~isempty(fig) && isvalid(fig)
    uiresume(fig);
    delete(fig);
end
end

function onProbabilityBarClick(~, ~)

if isempty(probAx) || ~isvalid(probAx)
    return;
end

cp = get(probAx, 'CurrentPoint');
clickedX = cp(1,1);

newIdx = round(clickedX);

if isnan(newIdx) || newIdx < 1 || newIdx > height(Tfull)
    return;
end

state.i = newIdx;

plotCurrent();
updateReviewNavigator();
drawnow;
end

function jumpToHighestProbUnreviewed()

if ~ismember('probSpike', Tfull.Properties.VariableNames)
    fprintf('No probSpike column found.\n');
    return;
end

idxUnrev = find(isnan(Tfull.reviewLabel));

if isempty(idxUnrev)
    fprintf('No unreviewed events remain.\n');
    return;
end

[~, imax] = max(Tfull.probSpike(idxUnrev));
state.i = idxUnrev(imax);

plotCurrent();
updateReviewNavigator();
drawnow;
end