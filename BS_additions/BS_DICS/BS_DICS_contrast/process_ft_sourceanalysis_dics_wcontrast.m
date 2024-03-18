function varargout = process_ft_sourceanalysis_dics_wcontrast(varargin )
% PROCESS_FT_SOURCEANALYSIS Call FieldTrip function ft_sourceanalysis (DICS)

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Vahab YoussofZadeh, Francois Tadel, 2021

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'FieldTrip: ft_sourceanalysis (wDICS) window contrast';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
%     sProcess.SubGroup    = 'Frequency';
sProcess.Index       = 357; %661
sProcess.Description = 'https://github.com/vyoussofzadeh/DICS-beamformer-for-Brainstorm';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data'};
sProcess.OutputTypes = {'data'};
sProcess.nInputs     = 2;
sProcess.nMinFiles   = 1;

% Option: Sensors selection
sProcess.options.sensortype.Comment = 'Sensor type:';
sProcess.options.sensortype.Type    = 'combobox_label';
sProcess.options.sensortype.Value   = {'MEG', {'MEG', 'MEG GRAD', 'MEG MAG', 'EEG', 'SEEG', 'ECOG'; ...
    'MEG', 'MEG GRAD', 'MEG MAG', 'EEG', 'SEEG', 'ECOG'}};
% Label: Time
sProcess.options.label1.Comment = '<BR><B>Time of interest:</B>';
sProcess.options.label1.Type    = 'label';
% Active time window
sProcess.options.poststim.Comment = 'Active (post-stim):';
sProcess.options.poststim.Type    = 'poststim';
sProcess.options.poststim.Value   = [];

% Label: Frequency
sProcess.options.label2.Comment = '<BR><B>Frequency of interset:</B>';
sProcess.options.label2.Type    = 'label';
% Enter the FOI in the data in Hz, eg, 22:
sProcess.options.foi.Comment = 'FOI:';
sProcess.options.foi.Type    = 'value';
sProcess.options.foi.Value   = {18, 'Hz', 0};
% Enter the Tapering frequency in the data in Hz, eg, 4:
sProcess.options.tpr.Comment = 'Tapering freq:';
sProcess.options.tpr.Type    = 'value';
sProcess.options.tpr.Value   = {4, 'Hz', 0};

% Label: Window
sProcess.options.label5.Comment = '<BR><B> Window of interset:</B>';
sProcess.options.label5.Type    = 'label';
sProcess.options.tlength.Comment = 'Length:';
sProcess.options.tlength.Type    = 'value';
sProcess.options.tlength.Value   = {0.3, 'ms',0};

sProcess.options.ovp.Comment = 'Overlap:';
sProcess.options.ovp.Type    = 'value';
sProcess.options.ovp.Value   =  {0.5, '0 to 1 (%100)', 1}; %{0.0,'0 to 1 (%200)',0.0};

sProcess.options.simaps.Comment = 'Save interval maps';
sProcess.options.simaps.Type    = 'checkbox';
sProcess.options.simaps.Value   = 1;

sProcess.options.avmaps.Comment = 'Save avg map';
sProcess.options.avmaps.Type    = 'checkbox';
sProcess.options.avmaps.Value   = 1;

% Label: Contrast
sProcess.options.label4.Comment = '<BR><B>Contrasting analysis:</B>'; % Contrast between pre and post, across trials
sProcess.options.label4.Type    = 'label';
% Contrast
sProcess.options.method.Comment = {'Subtraction (post-pre)', 'Permutation-stats (post-vs-pre)', 'Contrasting:'; ...
    'subtraction', 'permutation', ''};
sProcess.options.method.Type    = 'radio_linelabel';
sProcess.options.method.Value   = 'subtraction';
% ERDS
sProcess.options.erds.Comment = {'ERD', 'ERS', 'both', 'ERD/S effects:'; ...
    'erd', 'ers', 'both', ''};
sProcess.options.erds.Type    = 'radio_linelabel';
sProcess.options.erds.Value   = 'erd';
% Effects
sProcess.options.effect.Comment = {'abs', 'raw', 'Absolute/raw value of the contrast:'; ...
    'abs', 'raw', ''};
sProcess.options.effect.Type    = 'radio_linelabel';
sProcess.options.effect.Value   = 'abs';

% Label: TFR
sProcess.options.label3.Comment = '<BR><B>Time-freq response (used for data inspection):</B>';
sProcess.options.label3.Type    = 'label';
% Max frequency
sProcess.options.maxfreq.Comment = 'Max frequency (for TFR calculation):';
sProcess.options.maxfreq.Type    = 'value';
sProcess.options.maxfreq.Value   = {40, 'Hz', 0};
% Plot time-freq
sProcess.options.showtfr.Comment = 'Plot time-frequency figure for visual';
sProcess.options.showtfr.Type    = 'checkbox';
sProcess.options.showtfr.Value   = 1;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputsA, sInputsB) %#ok<DEFNU>
OutputFiles = {};
% Initialize FieldTrip
[isInstalled, errMsg] = bst_plugin('Install', 'fieldtrip');
if ~isInstalled
    bst_report('Error', sProcess, [], errMsg);
    return;
end

% ===== GET OPTIONS =====
% Inverse options
Method   = sProcess.options.method.Value;
Modality = sProcess.options.sensortype.Value{1};
ShowTfr  = sProcess.options.showtfr.Value;
MaxFreq  = sProcess.options.maxfreq.Value{1};
PostStim = sProcess.options.poststim.Value{1};
FOI      = sProcess.options.foi.Value{1};
TprFreq  = sProcess.options.tpr.Value{1};
TmpDir = bst_get('BrainstormTmpDir');
Overlap = sProcess.options.ovp.Value{1};
avmaps  = sProcess.options.avmaps.Value;
simaps  = sProcess.options.simaps.Value;
tlength = sProcess.options.tlength.Value{1};

% Load channel file
ChannelMat = in_bst_channel(sInputsA(1).ChannelFile);
% Get selected sensors
iChannels = channel_find(ChannelMat.Channel, Modality);
if isempty(iChannels)
    bst_report('Error', sProcess, sInputsA, ['Channels "' Modality '" not found in channel file.']);
    return;
end

% ===== LOAD: BAD CHANNELS =====
% Load bad channels from all the input files
isChannelGood = [];
for iInput = 1:length(sInputsA)
    DataFile = sInputsA(1).FileName;
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
% if isempty(iChannelsData)
%     bst_report('Error', sProcess, sInputs, 'All the selected channels are tagged as bad.');
%     return;
% elseif any(~isChannelGood)
%     bst_report('Info', sProcess, sInputs, ['Found ' num2str(length(find(~isChannelGood))) ' bad channels: ', sprintf('%s ', ChannelMat.Channel(~isChannelGood).Name)]);
% end

% ===== LOAD: HEADMODEL =====
% Get the study
sStudyChan = bst_get('ChannelFile', sInputsA(1).ChannelFile);
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
ftData = out_fieldtrip_data(sInputsA(1).FileName, ChannelMat, iChannelsData, 1);
ftData.trial = cell(1,length(sInputsA));
ftData.time = cell(1,length(sInputsA));
% Load all the trials

AllChannelFiles = unique({sInputsA.ChannelFile});
iChanInputs = find(ismember({sInputsA.ChannelFile}, AllChannelFiles{1}));
for iInput = 1:length(sInputsA)
    DataFile = sInputsA(iChanInputs(iInput)).FileName;
    DataMat = in_bst_data(DataFile);
    ftData.trial{iInput} = DataMat.F(iChannelsData,:);
    ftData.time{iInput} = DataMat.Time;
end
ftData_A = ftData;

% ===== LOAD: DATA =====
% Template FieldTrip structure for all trials
ftData = out_fieldtrip_data(sInputsB(1).FileName, ChannelMat, iChannelsData, 1);
ftData.trial = cell(1,length(sInputsB));
ftData.time = cell(1,length(sInputsB));
% Load all the trials

AllChannelFiles = unique({sInputsB.ChannelFile});
iChanInputs = find(ismember({sInputsB.ChannelFile}, AllChannelFiles{1}));
for iInput = 1:length(sInputsB)
    DataFile = sInputsB(iChanInputs(iInput)).FileName;
    DataMat = in_bst_data(DataFile);
    ftData.trial{iInput} = DataMat.F(iChannelsData,:);
    ftData.time{iInput} = DataMat.Time;
end
ftData_B = ftData;

%%
Overlap1 = 1-Overlap;
w1 = PostStim(1); l = tlength; ov = l.*Overlap1; j=1; wi=[];
while w1+l <= PostStim(2)
    wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
end
disp(wi)

%%
clear dics
for j=1:size(wi,1)
    
    disp(wi(j,:))
    % ===== FIELDTRIP: EPOCHING =====
    % Baseline
    cfg = [];
    % Post-stim
    cfg.toilim = wi(j,:); %PostStim;
    ep_data.pst = ft_redefinetrial(cfg, ftData_A);
    cfg.toilim = PostStim;
    ep_data.bsl = ft_redefinetrial(cfg, ftData_B);
    % Baseline + Post-stim
    cfg = [];
    ep_data.app = ft_appenddata(cfg, ep_data.bsl, ep_data.pst);
    
    %%
    switch sProcess.options.sensortype.Value{1}
        case {'EEG', 'SEEG', 'ECOG'}
            sens = ftData.elec;
        case {'MEG', 'MEG GRAD', 'MEG MAG'}
            sens = ftData.grad;
    end
    
    cfg = [];
    cfg.foilim = [FOI, FOI];
    if (FOI >= 4)
        cfg.taper = 'dpss';
        cfg.tapsmofrq = TprFreq;
    else
        cfg.taper = 'hanning';
        cfg.tapsmofrq = 1;
    end
    
    f_data.app = do_fft(cfg, ep_data.app); f_data.app.elec = sens;
    f_data.pst = do_fft(cfg, ep_data.pst); f_data.pst.elec = sens;
    f_data.bsl = do_fft(cfg, ep_data.bsl); f_data.bsl.elec = sens;
    
    %%
    % ===== SOURCE ANALYSIS =====
    switch Method
        case 'subtraction'
            cfg = [];
            cfg.method = 'dics';
            cfg.dics.lambda = '100%';
            cfg.sourcemodel  = ftLeadfield;
            cfg.frequency    = f_data.app.freq;
            cfg.headmodel = ftHeadmodel;
            cfg.dics.keepfilter = 'yes';
            cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
            sourceavg = ft_sourceanalysis(cfg, f_data.app);
            
            cfg = [];
            cfg.method = 'dics';
            cfg.dics.lambda = '0%';
            cfg.sourcemodel        = ftLeadfield;
            cfg.sourcemodel.filter = sourceavg.avg.filter;
            cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
            cfg.headmodel = ftHeadmodel;
            s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
            s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);
            
            switch sProcess.options.erds.Value
                case 'erd'
                    cfg = [];
                    cfg.parameter = 'pow';
                    cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
                    source_diff_dics = ft_math(cfg, s_data.pst, s_data.bsl);
                    source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
                    source_diff_dics.pow(source_diff_dics.pow>0)=0;
                case 'ers'
                    cfg = [];
                    cfg.parameter = 'pow';
                    cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
                    source_diff_dics = ft_math(cfg, s_data.pst, s_data.bsl);
                    source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
                    source_diff_dics.pow(source_diff_dics.pow<0)=0;
                case 'both'
                    cfg = [];
                    cfg.parameter = 'pow';
                    cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
                    source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
                    source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
            end
            
        case 'permutation'
            cfg = [];
            cfg.method = 'dics';
            cfg.dics.lambda = '100%';
            cfg.frequency    = f_data.app.freq;
            cfg.headmodel = ftHeadmodel;
            cfg.sourcemodel  = ftLeadfield;
            cfg.dics.keepfilter = 'yes';
            cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
            sourceavg = ft_sourceanalysis(cfg, f_data.app);
            
            cfg = [];
            cfg.method = 'dics';
            cfg.sourcemodel        = ftLeadfield;
            cfg.sourcemodel.filter = sourceavg.avg.filter;
            cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
            cfg.rawtrial = 'yes';
            cfg.headmodel = ftHeadmodel;
            s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
            s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);
            
            stat = do_source_stat_montcarlo(s_data);
            
            tmp = stat.stat;
            tmp2 = zeros(size(stat.pos,1),1);
            tmp2(stat.inside) = tmp;
            
            stats1  = stat;
            stats1.stat =  tmp2;
            stats1.mask = stat.inside;
            stats2 = stats1;
            stats2.stat(stats2.stat>0)=0;
            stats2.stat(isnan(stats2.stat))=0;
    end
    dics(j,:) = source_diff_dics.pow;
end

%%
DataMat.Time = wi(1,1):.01:wi(end,2);
n_verticies = size(dics,2);

% Initialize the output matrix with zeros
output_matrix = zeros(n_verticies,length(DataMat.Time));

% For each computed interval
for i = 1:j
    % Get the start and end indices
    start_idx = find(abs(DataMat.Time - wi(i,1)) < 1e-10);
    end_idx = find(abs(DataMat.Time - wi(i, 2)) < 1e-10);
    
    % Populate the output matrix with computed source activities
    % For this example, I'll just use random values, but you should use your actual computed values
    output_matrix(:, start_idx:end_idx) = dics(i,:)' + output_matrix(:, start_idx:end_idx);
end
D = output_matrix;

%%
% ===== SAVE RESULTS =====
% === CREATE OUTPUT STRUCTURE ===
bst_progress('text', 'Saving source file...');
bst_progress('inc', 1);
% Output study
if (length(sInputsA) == 1)
    iStudyOut = sInputsA(1).iStudy;
    RefDataFile = sInputs(iChanInputs(sInputsA)).FileName;
else
    [tmp, iStudyOut] = bst_process('GetOutputStudy', sProcess, sInputsA);
    RefDataFile = [];
end
% Create structure
ResultsMat = db_template('resultsmat');
ResultsMat.ImagingKernel = [];
switch Method
    case 'subtraction'
        switch sProcess.options.effect.Value
            case 'abs'
                source_diff_dics.pow = abs(D);
            case 'raw'
                source_diff_dics.pow = D;
        end
        ResultsMat.ImageGridAmp  = source_diff_dics.pow;
        ResultsMat.cfg           = source_diff_dics.cfg;
        
    case 'permutation'
        switch sProcess.options.effect.Value
            case 'abs'
                stats2.stat = abs((stats2.stat));
            case 'raw'
                stats2.stat = stats2.stat;
        end
        ResultsMat.ImageGridAmp  = stats2.stat;
        ResultsMat.cfg           = stat.cfg;
end
ResultsMat.nComponents   = 1;
ResultsMat.Function      = Method;
ResultsMat.Time          = DataMat.Time;
ResultsMat.DataFile      = RefDataFile;
ResultsMat.HeadModelFile = HeadModelFile;
ResultsMat.HeadModelType = HeadModelMat.HeadModelType;
ResultsMat.ChannelFlag   = DataMat.ChannelFlag;
ResultsMat.GoodChannel   = iChannelsData;
ResultsMat.SurfaceFile   = HeadModelMat.SurfaceFile;
ResultsMat.nAvg          = DataMat.nAvg;
ResultsMat.Leff          = DataMat.Leff;

if isfield(sProcess.options, 'comment')
    ResultsMat.Comment = sProcess.options.comment.Value;
else
    ResultsMat.Comment       = ['wDICS: ' Method, ' ',num2str(FOI),'Hz ', sprintf('%1.3fs-%1.3fs', PostStim), ' WinL ', num2str(tlength), 's OverL ', num2str(100.*Overlap), '%'];
end
switch lower(ResultsMat.HeadModelType)
    case 'volume'
        ResultsMat.GridLoc    = HeadModelMat.GridLoc;
    case 'surface'
        ResultsMat.GridLoc    = [];
    case 'mixed'
        ResultsMat.GridLoc    = HeadModelMat.GridLoc;
        ResultsMat.GridOrient = HeadModelMat.GridOrient;
end
ResultsMat = bst_history('add', ResultsMat, 'compute', ['ft_sourceanalysis: ' Method ' ' Modality ' ']);

% === SAVE OUTPUT FILE ===
% Output filename
OutputDir = bst_fileparts(file_fullpath(DataFile));
ResultFile = bst_process('GetNewFilename', OutputDir, ['results_', Method, '_', Modality]);
% Save new file structure
bst_save(ResultFile, ResultsMat, 'v6');

% ===== REGISTER NEW FILE =====
% Create new results structure
newResult = db_template('results');
newResult.Comment       = ResultsMat.Comment;
newResult.FileName      = file_short(ResultFile);
newResult.DataFile      = ResultsMat.DataFile;
newResult.isLink        = 0;
newResult.HeadModelType = ResultsMat.HeadModelType;
% Get output study
sStudyOut = bst_get('Study', iStudyOut);
% Add new entry to the database
iResult = length(sStudyOut.Result) + 1;
sStudyOut.Result(iResult) = newResult;
% Update Brainstorm database
bst_set('Study', iStudyOut, sStudyOut);
% Store output filename
OutputFiles{end+1} = newResult.FileName;
% Expand data node
panel_protocols('SelectNode', [], newResult.FileName);

% Save database
db_save();
% Hide progress bar
bst_progress('stop');

end

%% ===== TIME-FREQUENCY =====
function [time_of_interest,freq_of_interest] = do_tfr_plot(cfg_main, tfr)
% First compute the average over trials:
cfg = [];
freq_avg = ft_freqdescriptives(cfg, tfr);

% And baseline-correct the average:
cfg = [];
cfg.baseline = cfg_main.baselinetime;
cfg.baselinetype = 'db'; % Use decibel contrast here
freq_avg_bsl = ft_freqbaseline(cfg, freq_avg);

freq_avg_bsl.powspctrm(isnan(freq_avg_bsl.powspctrm))=0;
meanpow = squeeze(mean(freq_avg_bsl.powspctrm, 1));

tim_interp = linspace(cfg_main.toi(1), cfg_main.toi(2), 512);
freq_interp = linspace(1, cfg_main.fmax, 512);

% We need to make a full time/frequency grid of both the original and
% interpolated coordinates. Matlab's meshgrid() does this for us:
[tim_grid_orig, freq_grid_orig] = meshgrid(tfr.time, tfr.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

% And interpolate:
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow, tim_grid_interp, freq_grid_interp, 'spline');

% while n==1
pow_interp1  = pow_interp(50:end,50:end);
tim_interp1  = tim_interp(50:end);
freq_interp1 = freq_interp(50:end);


[~,idx] = min(pow_interp1(:));
[row,col] = ind2sub(size(pow_interp1),idx);

time_of_interest = tim_interp1(col);
freq_of_interest = freq_interp1(row);

timind = nearest(tim_interp, time_of_interest);
freqind = nearest(freq_interp, freq_of_interest);
pow_at_toi = pow_interp(:,timind);
pow_at_foi = pow_interp(freqind,:);

% Plot figure
if cfg_main.plotflag
    figure();
    ax_main  = axes('Position', [0.1 0.2 0.55 0.55]);
    ax_right = axes('Position', [0.7 0.2 0.1 0.55]);
    ax_top   = axes('Position', [0.1 0.8 0.55 0.1]);
    
    axes(ax_main);
    imagesc(tim_interp, freq_interp, pow_interp);
    % note we're storing a handle to the image im_main, needed later on
    xlim([cfg_main.toi(1), cfg_main.toi(2)]);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    clim = max(abs(meanpow(:)));
    caxis([-clim clim]);
    % colormap(brewermap(256, '*RdYlBu'));
    hold on;
    plot(zeros(size(freq_interp)), freq_interp, 'k:');
    
    axes(ax_top);
    area(tim_interp, pow_at_foi, ...
        'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
    xlim([cfg_main.toi(1), cfg_main.toi(2)]);
    ylim([-clim clim]);
    box off;
    ax_top.XTickLabel = [];
    ylabel('Power (dB)');
    hold on;
    plot([0 0], [-clim clim], 'k:');
    
    axes(ax_right);
    area(freq_interp, pow_at_toi,...
        'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
    view([270 90]); % this rotates the plot
    ax_right.YDir = 'reverse';
    ylim([-clim clim]);
    box off;
    ax_right.XTickLabel = [];
    ylabel('Power (dB)');
    
    h = colorbar(ax_main, 'manual', 'Position', [0.85 0.2 0.05 0.55]);
    ylabel(h, 'Power vs baseline (dB)');
    
    % Main plot:
    axes(ax_main);
    plot(ones(size(freq_interp))*time_of_interest, freq_interp,...
        'Color', [0 0 0 0.1], 'LineWidth', 3);
    plot(tim_interp, ones(size(tim_interp))*freq_of_interest,...
        'Color', [0 0 0 0.1], 'LineWidth', 3);
    
    % Marginals:
    axes(ax_top);
    plot([time_of_interest time_of_interest], [0 clim],...
        'Color', [0 0 0 0.1], 'LineWidth', 3);
    axes(ax_right);
    hold on;
    plot([freq_of_interest freq_of_interest], [0 clim],...
        'Color', [0 0 0 0.1], 'LineWidth', 3);
end
end

function [freq,ff, psd,tapsmofrq] = do_fft(cfg_mian, data)
cfg              = [];
cfg.method       = 'mtmfft';
cfg.output       = 'fourier';
cfg.keeptrials   = 'yes';
cfg.foilim       = cfg_mian.foilim;
cfg.tapsmofrq    = cfg_mian.tapsmofrq;
cfg.taper        = cfg_mian.taper;
cfg.pad          = 4;
freq             = ft_freqanalysis(cfg, data);
psd = squeeze(mean(mean(abs(freq.fourierspctrm),2),1));
ff = linspace(1, cfg.foilim(2), length(psd));

tapsmofrq = cfg.tapsmofrq;
end

function stat = do_source_stat_montcarlo(s_data)
cfg = [];
cfg.parameter        = 'pow';
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'fdr';
cfg.clusteralpha     = 0.001;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 5000;

ntrials                       = numel(s_data.bsl.trial);
design                        = zeros(2,2*ntrials);
design(1,1:ntrials)           = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials)           = 1:ntrials;
design(2,ntrials+1:2*ntrials) = 1:ntrials;

cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;
stat         = ft_sourcestatistics(cfg,s_data.pst,s_data.bsl);
end


function [data] = ft_redefinetrial(cfg, data)

% FT_REDEFINETRIAL allows you to adjust the time axis of your data, i.e. to
% change from stimulus-locked to response-locked. Furthermore, it allows
% you to select a time window of interest, or to resegment your long trials
% into shorter fragments.
%
% Use as
%   [data] = ft_redefinetrial(cfg, data)
% where the input data should correspond to the output of FT_PREPROCESSING and the
% configuration should be specified as explained below. Note that some options are
% mutually exclusive. If you want to use both,  you neew two calls to this function
% to avoid confusion about the order in which they are applied.
%
% For selecting a subset of trials you can specify
%   cfg.trials    = 'all' or a selection given as a 1xN vector (default = 'all')
%
% For selecting trials with a minimum length you can specify
%   cfg.minlength = length in seconds, can be 'maxperlen' (default = [])
%
% For realiging the time axes of all trials to a new reference time
% point (i.e. change the definition for t=0) you can use the following
% configuration option
%   cfg.offset    = single number or Nx1 vector, expressed in samples relative to current t=0
%
% For selecting a specific subsection of (i.e. cut out a time window
% of interest) you can select a time window in seconds that is common
% in all trials
%   cfg.toilim    = [tmin tmax] to specify a latency window in seconds, can be Nx2 vector
%
% Alternatively you can specify the begin and end sample in each trial
%   cfg.begsample = single number or Nx1 vector, expressed in samples relative to the start of the input trial
%   cfg.endsample = single number or Nx1 vector, expressed in samples relative to the start of the input trial
%
% Alternatively you can specify a new trial definition, expressed in
% samples relative to the original recording
%   cfg.trl       = Nx3 matrix with the trial definition, see FT_DEFINETRIAL
%
% Alternatively you can specify the data to be cut into (non-)overlapping
% segments, starting from the beginning of each trial. This may lead to loss
% of data at the end of the trials
%   cfg.length    = number (in seconds) that specifies the length of the required snippets
%   cfg.overlap   = number between 0 and 1 (exclusive) specifying the fraction of overlap between snippets (0 = no overlap)
%
% Alternatively you can merge or stitch pseudo-continuous segmented data back into a
% continuous representation. This requires that the data has a valid sampleinfo field
% and that there are no jumps in the signal in subsequent trials (e.g. due to
% filtering or demeaning). If there are missing segments (e.g. due to artifact
% rejection), the output data will have one trial for each section where the data is
% continuous.
%   cfg.continuous = 'yes'
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_DEFINETRIAL, FT_RECODEEVENT, FT_PREPROCESSING

% Copyright (C) 2006-2021, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
    return
end

% store original datatype
dtype = ft_datatype(data);

% deal with the special case of timelock rpt_chan_time with 1 trial
oneRptTimelock = (strcmp(dtype, 'timelock') && ...
    strcmp(data.dimord, 'rpt_chan_time') && ...
    size(data.trial, 1) == 1);

% check if the input data is valid for this function, this will convert it to raw if needed
data = ft_checkdata(data, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'trial'}); % prevent accidental typos, see issue 1729

% set the defaults
cfg.offset       = ft_getopt(cfg, 'offset',     []);
cfg.toilim       = ft_getopt(cfg, 'toilim',     []);
cfg.begsample    = ft_getopt(cfg, 'begsample',  []);
cfg.endsample    = ft_getopt(cfg, 'endsample',  []);
cfg.minlength    = ft_getopt(cfg, 'minlength',  []);
cfg.trials       = ft_getopt(cfg, 'trials',     'all', 1);
cfg.feedback     = ft_getopt(cfg, 'feedback',   'yes');
cfg.trl          = ft_getopt(cfg, 'trl',        []);
cfg.length       = ft_getopt(cfg, 'length',     []);
cfg.overlap      = ft_getopt(cfg, 'overlap',    0);
cfg.continuous   = ft_getopt(cfg, 'continuous', 'no');

% select trials of interest
if ~strcmp(cfg.trials, 'all')
    if islogical(cfg.trials)
        ft_info('selecting %d trials\n', sum(cfg.trials));
    else
        ft_info('selecting %d trials\n', length(cfg.trials));
    end
    
    % select trials of interest
    tmpcfg = keepfields(cfg, {'trials', 'showcallinfo', 'trackcallinfo', 'trackconfig', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
    data   = ft_selectdata(tmpcfg, data);
    % restore the provenance information
    [cfg, data] = rollback_provenance(cfg, data);
    
    if length(cfg.offset)>1 && length(cfg.offset)~=length(cfg.trials)
        cfg.offset = cfg.offset(cfg.trials);
    end
    if length(cfg.begsample)>1 && length(cfg.begsample)~=length(cfg.trials)
        cfg.begsample = cfg.begsample(cfg.trials);
    end
    if length(cfg.endsample)>1 && length(cfg.endsample)~=length(cfg.trials)
        cfg.endsample = cfg.endsample(cfg.trials);
    end
end

Ntrial = numel(data.trial);

% check the input arguments, only one method for processing is allowed
numoptions = ~isempty(cfg.toilim) + ~isempty(cfg.offset) + (~isempty(cfg.begsample) || ~isempty(cfg.endsample)) + ~isempty(cfg.trl) + ~isempty(cfg.length) + istrue(cfg.continuous);
if numoptions>1
    ft_error('you should specify only one of the options for redefining the data segments');
end
if numoptions==0 && isempty(cfg.minlength) && strcmp(cfg.trials, 'all')
    ft_error('you should specify at least one configuration option');
end

% start processing
if ~isempty(cfg.toilim)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select a latency window from each trial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if numel(cfg.toilim) == 2
        % specified as single [tstart tend] vector
        % expand into Ntrial X 2
        cfg.toilim = repmat(cfg.toilim(:)', Ntrial, 1);
    end
    
    begsample = zeros(Ntrial,1);
    endsample = zeros(Ntrial,1);
    skiptrial = false(Ntrial,1);
    for i=1:Ntrial
        if cfg.toilim(i,1)>data.time{i}(end) || cfg.toilim(i,2)<data.time{i}(1)
            begsample(i) = nan;
            endsample(i) = nan;
            skiptrial(i) = true;
        else
            begsample(i) = nearest(data.time{i}, cfg.toilim(i,1));
            endsample(i) = nearest(data.time{i}, cfg.toilim(i,2));
            data.trial{i} = data.trial{i}(:, begsample(i):endsample(i));
            data.time{i}  = data.time{i} (   begsample(i):endsample(i));
        end
    end
    
    % also correct the sample information
    if isfield(data, 'sampleinfo')
        data.sampleinfo(:, 1) = data.sampleinfo(:, 1) + begsample - 1;
        data.sampleinfo(:, 2) = data.sampleinfo(:, 1) + endsample - begsample;
    end
    
    data.time     = data.time(~skiptrial);
    data.trial    = data.trial(~skiptrial);
    if isfield(data, 'sampleinfo'),  data.sampleinfo  = data.sampleinfo(~skiptrial, :); end
    if isfield(data, 'trialinfo'),   data.trialinfo   = data.trialinfo(~skiptrial, :);  end
    ft_info('removing %d trials in which no data was selected\n', sum(skiptrial));
    
elseif ~isempty(cfg.offset)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % shift the time axis from each trial
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    offset = cfg.offset(:);
    offset = round(offset); % this is in samples and hence it must be expressed as integers
    if length(cfg.offset)==1
        offset = repmat(offset, Ntrial, 1);
    end
    for i=1:Ntrial
        data.time{i} = data.time{i} + offset(i)/data.fsample;
    end
    
elseif ~isempty(cfg.begsample) || ~isempty(cfg.endsample)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select a latency window from each trial based on begin and/or end sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begsample = cfg.begsample(:);
    endsample = cfg.endsample(:);
    if length(begsample)==1
        begsample = repmat(begsample, Ntrial, 1);
    end
    if length(endsample)==1
        endsample = repmat(endsample, Ntrial, 1);
    end
    for i=1:Ntrial
        data.trial{i} = data.trial{i}(:, begsample(i):endsample(i));
        data.time{i}  = data.time{i} (   begsample(i):endsample(i));
    end
    
    % also correct the sampleinfo
    if isfield(data, 'sampleinfo')
        sampleinfo = data.sampleinfo(:, 1);
        data.sampleinfo(:, 1) = sampleinfo(:, 1) + begsample - 1;
        data.sampleinfo(:, 2) = sampleinfo(:, 1) + endsample - 1;
    end
    
elseif ~isempty(cfg.trl)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select new trials from the existing data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ischar(cfg.trl)
        % load the trial information from file
        newtrl = loadvar(cfg.trl, 'trl');
    else
        newtrl = cfg.trl;
    end
    
    % ensure that sampleinfo is present, otherwise ft_fetch_data will crash
    data = ft_checkdata(data, 'hassampleinfo', 'yes');
    
    % make a copy of the old data
    dataold = data;
    
    % make the header
    hdr = ft_fetch_header(dataold);
    
    % start with a completely new data structure
    data          = keepfields(dataold, {'cfg' 'fsample' 'label' 'topo' 'topolabel' 'unmixing' 'mixing' 'grad' 'elec' 'opto'}); % account for all potential fields to be copied over
    data.hdr      = hdr;
    data.trial    = cell(1,size(newtrl,1));
    data.time     = cell(1,size(newtrl,1));
    
    if isfield(dataold, 'trialinfo')
        ft_warning('Original data has trialinfo, using user-specified trialinfo instead');
    end
    
    if ~istable(newtrl)
        begsample = newtrl(:,1);
        endsample = newtrl(:,2);
        offset    = newtrl(:,3);
    else
        begsample = newtrl.begsample;
        endsample = newtrl.endsample;
        offset    = newtrl.offset;
    end
    trllength = endsample - begsample + 1;
    
    for iTrl=1:size(newtrl, 1)
        
        data.trial{iTrl} = ft_fetch_data(dataold, 'header', hdr, 'begsample', begsample(iTrl), 'endsample', endsample(iTrl), 'chanindx', 1:hdr.nChans, 'skipcheckdata', 1);
        data.time{iTrl}  = offset2time(offset(iTrl), dataold.fsample, trllength(iTrl));
        
        % The following ensures correct handling of trialinfo.
        
        % Determine which old trials are present in new trials
        iTrlorig = find(begsample(iTrl) <= dataold.sampleinfo(:,2) & endsample(iTrl) >= dataold.sampleinfo(:,1));
        
        if size(newtrl,2)>3 % In case user specified additional trialinfo
            data.trialinfo(iTrl,:) = newtrl(iTrl,4:end);
        elseif isfield(dataold,'trialinfo') % If old data has trialinfo
            if (numel(iTrlorig) == 1 ...      % only 1 old trial to copy trialinfo from, or
                    || size(unique(dataold.trialinfo(iTrlorig,:),'rows'),1)) ... % all old trialinfo rows are identical
                    && ~any(diff(dataold.sampleinfo(:,1))<=0) % and the trials are consecutive segments
                data.trialinfo(iTrl,:) = dataold.trialinfo(iTrlorig(1),:);
            else
                ft_error('Old trialinfo cannot be combined into new trialinfo, please specify trialinfo in cfg.trl(:,4)');
            end
        end
    end % for iTrl
    
    % adjust the sampleinfo in the output
    if isfield(dataold, 'sampleinfo')
        % adjust the sample information
        data.sampleinfo  = [begsample endsample];
    end
    
elseif ~isempty(cfg.length)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cut the existing trials into segments of the specified length
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data = ft_checkdata(data, 'hassampleinfo', 'yes');
    
    % create dummy trl-matrix and recursively call ft_redefinetrial
    nsmp    = round(cfg.length*data.fsample);
    nshift  = round((1-cfg.overlap)*nsmp);
    
    newtrl = zeros(0,4);
    for k = 1:numel(data.trial)
        begsample = data.sampleinfo(k,1);
        endsample = data.sampleinfo(k,2);
        offset    = time2offset(data.time{k}, data.fsample);
        thistrl   = (begsample:nshift:(endsample+1-nsmp))';
        if ~isempty(thistrl) % the trial might be too short
            thistrl(:,2) = thistrl(:,1) + nsmp - 1;
            thistrl(:,3) = thistrl(:,1) + offset - thistrl(1,1);
            thistrl(:,4) = k; % keep the trial number in the 4th column, this is needed further down
            newtrl = cat(1, newtrl, thistrl);
        end
    end
    clear begsample endsample offset
    
    tmpcfg = keepfields(cfg, {'feedback', 'showcallinfo', 'trackcallinfo', 'trackconfig', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
    tmpcfg.trl = newtrl;
    
    if isfield(data, 'trialinfo') && ~istable(data.trialinfo)
        % replace the trial number with the original trial information
        tmpcfg.trl = [newtrl(:,1:3) data.trialinfo(newtrl(:,4),:)];
    elseif isfield(data, 'trialinfo') && istable(data.trialinfo)
        % construct the trl matrix as a table
        begsample = newtrl(:,1);
        endsample = newtrl(:,2);
        offset    = newtrl(:,3);
        tmpcfg.trl = [table(begsample, endsample, offset) data.trialinfo(newtrl(:,4),:)];
    elseif ~isfield(data, 'trialinfo')
        % discard the trial number
        tmpcfg.trl = newtrl(:,1:3);
    end
    
    data   = removefields(data, {'trialinfo'}); % these are in the additional columns of tmpcfg.trl
    data   = ft_redefinetrial(tmpcfg, data);
    % restore the provenance information
    [cfg, data] = rollback_provenance(cfg, data);
    
elseif istrue(cfg.continuous)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify consecutive segments that can be glued back together
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data = ft_checkdata(data, 'hassampleinfo', 'yes');
    
    boolvec = artifact2boolvec(data.sampleinfo);
    newtrl = boolvec2trl(boolvec);
    
    % In general: An offset of 0 means that the first sample of the trial corresponds
    % to the trigger. A positive offset indicates that the first sample is later than
    % the trigger.
    
    % here we want to use the start of the recording as t=0
    newtrl(:,3) = newtrl(:,1) - 1;
    
    tmpcfg = keepfields(cfg, {'feedback', 'showcallinfo', 'trackcallinfo', 'trackconfig', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
    tmpcfg.trl = newtrl;
    
    data   = removefields(data, {'trialinfo'}); % the trialinfo does not apply any more
    data   = ft_redefinetrial(tmpcfg, data);
    % restore the provenance information
    [cfg, data] = rollback_provenance(cfg, data);
    
end % processing the realignment or data selection

if ~isempty(cfg.minlength)
    Ntrial    = length(data.trial);
    trllength = zeros(Ntrial, 1);
    % determine the length of each trial
    for i=1:Ntrial
        trllength(i) = size(data.trial{i},2) * 1/data.fsample; % this the the DURATION of the selected samples
    end
    if ischar(cfg.minlength) && strcmp(cfg.minlength, 'maxperlen')
        minlength = max(trllength);
    else
        minlength = cfg.minlength;
    end
    % remove trials that are too short
    skiptrial = (trllength<minlength);
    %if ~isempty(trl), trl = trl(~skiptrial,:); end
    data.time  = data.time(~skiptrial);
    data.trial = data.trial(~skiptrial);
    if isfield(data, 'sampleinfo'), data.sampleinfo  = data.sampleinfo(~skiptrial, :); end
    if isfield(data, 'trialinfo'),  data.trialinfo   = data.trialinfo (~skiptrial, :); end
    ft_info('removing %d trials that are too short\n', sum(skiptrial));
end

% convert back to input type if necessary
switch dtype
    case 'timelock'
        data = ft_checkdata(data, 'datatype', 'timelock');
        if oneRptTimelock
            % deal with the special case of rpt_chan_time timelock data with one
            % repetition
            data.trial = reshape(data.avg, [1 size(data.avg)]);
            data.dimord = 'rpt_chan_time';
        end
    otherwise
        % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
end
