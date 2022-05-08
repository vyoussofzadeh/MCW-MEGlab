function varargout = process_ft_reject_trials(varargin )
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
% Authors: Vahab YoussofZadeh, 2022, update, 02/15/2022

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'FieldTrip: process_ft_reject_trials';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Pre-process';
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

% Rejection method
sProcess.options.label2.Comment = '<BR><B>Rejection method:</B>';
sProcess.options.label2.Type    = 'label';

sProcess.options.mrej.Comment = {'Automatic', 'Manual', 'Rejection method:'; ...
    'Auto', 'Manual', ''};
sProcess.options.mrej.Type    = 'radio_linelabel';
sProcess.options.mrej.Value   = 'RejMethod';

sProcess.options.arej.Comment = 'Automatic rejection';
sProcess.options.arej.Type    = 'value';
sProcess.options.arej.Value   = {0.9, 'threshold', 1};

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
OutputFiles = {};
% Initialize FieldTrip
[isInstalled, errMsg] = bst_plugin('Install', 'fieldtrip');
if ~isInstalled
    bst_report('Error', sProcess, [], errMsg);
    return;
end

% ===== GET OPTIONS =====
Modality = sProcess.options.sensortype.Value{1};
arej = sProcess.options.arej.Value;
mrej = sProcess.options.mrej.Value;

% Output folder (for figures and FieldTrip structures)
TmpDir = bst_get('BrainstormTmpDir');
% Progress bar
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
[~, ~, iChannelsData] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);

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

%%
ftData1 = ftData;
ftData1.grad.chantype = ftData1.grad.chantype';

%%
close all;
switch mrej
    case 'Manual'
        ftData2 = ftData1; ftData2.avg = [];       
        cfg = [];
        cfg.metric = 'kurtosis';  % use by default kurtosis method
        cfg.latency = [ftData1.time{1}(1),ftData1.time{1}(end)];
        r_data   = ft_rejectvisual(cfg, ftData2);
        
        % Bad trial info
        data = ft_checkdata(ftData2, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');
        btrlsample = r_data.cfg.artfctdef.summary.artifact;
        for l=1:size(btrlsample,1)
            btrl(l,:) = find(ismember(data.sampleinfo(:,1),btrlsample(l))==1);
        end
        disp('BAD TRIAL:')
        disp(btrl)
        
    case 'Auto'
        cfg = [];
        cfg.pflag = 1; % Yes:1, No:2
        cfg.saveflag = 2; % Yes:1, No:2
        cfg.savepath = [];
        cfg.latency = [ftData1.time{1}(1),ftData1.time{1}(end)];
        cfg.rejectpercentage = arej{1};
        cfg.rbadtrl = 1;
        cfg.rbadsen = 0;
        [r_data, report] = do_artifactreject(cfg, ftData1);
        if cfg.rbadtrl == 1
            disp('BAD TRIAL:')
            disp(report.btrl')
        end
end
%%

%%
disp('np = 0, yes = 1');
disp('apply the correction):')
bic = input(['']);
close all;

if bic == 1   
    %-
    ProtocolInfo = bst_get('ProtocolInfo');
    [datapath,~] = fileparts(DataFile);
    [~,srcPath] = fileparts(datapath);
    srcDir  = bst_fullfile(ProtocolInfo.STUDIES, [srcPath,'_arj']);
    
    cd(fullfile(ProtocolInfo.STUDIES, datapath))
    for iInput = 1:length(sInputs)
        DataFile = sInputs(iInput).FileName;
        D = load(file_fullpath(DataFile));
        [a,FileName] = fileparts(DataFile);
        FileName_new = [FileName, '_arj'];
        tkz = tokenize(D.Comment, ' ');
        D.Comment = [tkz{1} '_trl ',tkz{2}];
        D.F(iChannelsData,:) = r_data.trial{iInput};
        save(FileName_new, '-struct', 'D');        
    end
    disp('done, reload datafile!')
    
else
    disp('no correction was done')
end
db_save();
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

function [r_data, report] = do_artifactreject(cfg_main, dat)

% disp('Identifying bad channels/sensors and bad trials ...');
% if exist(cfg_main.savepath, 'file') == 2
%     load(cfg_main.savepath)
% else

%% kurtosis
cfg = [];
cfg.trials = 'all';
cfg.metric = 'kurtosis';
cfg.channel = 'all';
cfg.latency = cfg_main.latency;
[level,info] = do_compute_metric(cfg,dat);
%     metric.kurt = level;
info.pflag = cfg_main.pflag;
[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = do_plot_chantrl(info,level);

thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel

%% zvalue
cfg = [];
cfg.trials = 'all';
cfg.metric = 'zvalue';
cfg.channel = 'all';
cfg.latency = cfg_main.latency;
[level,info] = do_compute_metric(cfg,dat);
%     metric.kurt = level;
info.pflag = cfg_main.pflag;
[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = do_plot_chantrl(info,level);

thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel

%% Var
cfg = [];
cfg.trials = 'all';
cfg.metric = 'var';
cfg.channel = 'all';
cfg.latency = cfg_main.latency;
[level,info] = do_compute_metric(cfg,dat);

info.pflag = cfg_main.pflag;
[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = do_plot_chantrl(info,level);

thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel

%%
%     btrl_all = unique([btrl.kurtosis,btrl.var,btrl.zvalue]);
bch_all = unique([bch.kurtosis;bch.var;bch.zvalue]);

disp('Bad channels:')
for i=1:length(bch_all)
    bch_all_label_disp{i,:} = dat.label{bch_all(i)};
    bch_all_label{i,:} = ['-',dat.label{bch_all(i)}];
end
disp(bch_all_label_disp);

%% Removing bad trials
if cfg_main.rbadtrl == 1
    btrl_all = unique([btrl.kurtosis,btrl.var,btrl.zvalue]);
    disp('Bad trials:')
    disp(btrl_all);
    cfg = [];
    cfg.trials = find(~ismember(1:length(dat.trial),btrl_all));
    dat = ft_selectdata(cfg, dat);
    report.btrl = btrl_all;
end

%% Removing bad sensors
if cfg_main.rbadsen == 1
    if length(bch_all) < 10
        cfg = [];
        cfg.channel = ['all';bch_all_label];
        dat = ft_selectdata(cfg, dat);
    else
        warning('too many bad sensors w> 10, rejection skipped');
    end
    report.bchan = bch_all_label_disp;
end

r_data = dat;
if cfg_main.saveflag ==1
    %         save(cfg_main.savepath, 'r_data', 'report','-v7.3');
end
end

function [level, info] = do_compute_metric(cfg, data)
% SUBFUNCTION for ft_rejectvisual

% determine the initial selection of trials
ntrl = length(data.trial);
if isequal(cfg.trials, 'all') % support specification like 'all'
    cfg.trials = 1:ntrl;
end
trlsel = false(1, ntrl);
trlsel(cfg.trials) = true;

% determine the initial selection of channels
nchan = length(data.label);
cfg.channel = ft_channelselection(cfg.channel, data.label); % support specification like 'all'
chansel = false(1, nchan);
chansel(match_str(data.label, cfg.channel)) = true;

% compute the sampling frequency from the first two timepoints
fsample = 1/mean(diff(data.time{1}));

% select the specified latency window from the data
% here it is done BEFORE filtering and metric computation
for i=1:ntrl
    begsample = nearest(data.time{i}, cfg.latency(1));
    endsample = nearest(data.time{i}, cfg.latency(2));
    data.time{i} = data.time{i}(begsample:endsample);
    data.trial{i} = data.trial{i}(:, begsample:endsample);
end

% compute the offset from the time axes
offset = zeros(ntrl, 1);
for i=1:ntrl
    offset(i) = do_time2offset(data.time{i}, fsample);
end

% set up guidata info
info                = [];
info.data           = data;
info.cfg            = cfg;
info.metric         = cfg.metric;
info.previousmetric = 'none';
info.level          = nan(nchan, ntrl);
info.ntrl           = ntrl;
info.nchan          = nchan;
info.trlsel         = trlsel;
info.chansel        = chansel;
info.fsample        = fsample;
info.offset         = offset;
info.quit           = 0;


% update_log(cfg.output_box, 'Computing metric...');
% ft_progress('init', cfg.cfg.feedback, 'computing metric');
level = zeros(info.nchan, info.ntrl);
if strcmp(info.metric, 'zvalue') || strcmp(info.metric, 'maxzvalue')
    % cellmean and cellstd (see ft_denoise_pca) would work instead of for-loops, but they are too memory-intensive
    runsum = zeros(info.nchan, 1);
    runss  = zeros(info.nchan, 1);
    runnum = 0;
    for i=1:info.ntrl
        dat = do_preproc(info.data.trial{i}, info.data.label, do_offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)),[]); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
        runsum = runsum + nansum(dat, 2);
        runss  = runss  + nansum(dat.^2, 2);
        runnum = runnum + sum(isfinite(dat), 2);
    end
    mval = runsum./runnum;
    sd   = sqrt(runss./runnum - (runsum./runnum).^2);
end
for i=1:info.ntrl
    ft_progress(i/info.ntrl, 'computing metric %d of %d\n', i, info.ntrl);
    dat = do_preproc(info.data.trial{i}, info.data.label, do_offset2time(info.offset(i), info.fsample, size(info.data.trial{i}, 2)),[]); % not entirely sure whether info.data.time{i} is correct, so making it on the fly
    switch info.metric
        case 'var'
            level(:, i) = nanstd(dat, [], 2).^2;
        case 'min'
            level(:, i) = nanmin(dat, [], 2);
        case 'max'
            level(:, i) = nanmax(dat, [], 2);
        case 'maxabs'
            level(:, i) = nanmax(abs(dat), [], 2);
        case 'range'
            level(:, i) = nanmax(dat, [], 2) - nanmin(dat, [], 2);
        case 'kurtosis'
            level(:, i) = kurtosis(dat, [], 2);
        case '1/var'
            level(:, i) = 1./(nanstd(dat, [], 2).^2);
        case 'zvalue'
            level(:, i) = nanmean( (dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , 2);
        case 'maxzvalue'
            level(:, i) = nanmax( ( dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , [], 2);
        otherwise
            ft_error('unsupported method');
    end
end
end

function offset = do_time2offset(ttime, fsample)

% TIME2OFFSET converts a time-axis of a trial into the offset in samples
% according to the definition from DEFINETRIAL
%
% Use as
%   [offset] = time2offset(time, fsample)
%
% The trialdefinition "trl" is an Nx3 matrix. The first column contains
% the sample-indices of the begin of the trial relative to the begin
% of the raw data , the second column contains the sample_indices of
% the end of the trials, and the third column contains the offset of
% the trigger with respect to the trial. An offset of 0 means that
% the first sample of the trial corresponds to the trigger. A positive
% offset indicates that the first sample is later than the triger, a
% negative offset indicates a trial beginning before the trigger.

% Copyright (C) 2005, Robert Oostenveld
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

offset = round(ttime(1)*fsample);
end

function ttime = do_offset2time(offset, fsample, nsamples)

% OFFSET2TIME converts the offset of a trial definition into a time-axis
% according to the definition from DEFINETRIAL
%
% Use as
%   [time] = offset2time(offset, fsample, nsamples)
%
% The trialdefinition "trl" is an Nx3 matrix. The first column contains
% the sample-indices of the begin of the trial relative to the begin
% of the raw data , the second column contains the sample_indices of
% the end of the trials, and the third column contains the offset of
% the trigger with respect to the trial. An offset of 0 means that
% the first sample of the trial corresponds to the trigger. A positive
% offset indicates that the first sample is later than the triger, a
% negative offset indicates a trial beginning before the trigger.

% Copyright (C) 2005, Robert Oostenveld
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

% ensure that these are not integers
offset   = double(offset);
nsamples = double(nsamples);

ttime = (offset + (0:(nsamples-1)))/fsample;
end

function [dat, label, ttime, cfg] = do_preproc(dat, label, ttime, cfg, begpadding, endpadding)

% PREPROC applies various preprocessing steps on a piece of EEG/MEG data
% that already has been read from a data file.
%
% This function can serve as a subfunction for all FieldTrip modules that
% want to preprocess the data, such as PREPROCESSING, ARTIFACT_XXX,
% TIMELOCKANALYSIS, etc. It ensures consistent handling of both MEG and EEG
% data and consistency in the use of all preprocessing configuration
% options.
%
% Use as
%   [dat, label, time, cfg] = preproc(dat, label, time, cfg, begpadding, endpadding)
%
% The required input arguments are
%   dat         Nchan x Ntime data matrix
%   label       Nchan x 1 cell-array with channel labels
%   time        Ntime x 1 vector with the latency in seconds
%   cfg         configuration structure, see below
% and the optional input arguments are
%   begpadding  number of samples that was used for padding (see below)
%   endpadding  number of samples that was used for padding (see below)
%
% The output is
%   dat         Nchan x Ntime data matrix
%   label       Nchan x 1 cell-array with channel labels
%   time        Ntime x 1 vector with the latency in seconds
%   cfg         configuration structure, optionally with extra defaults set
%
% Note that the number of input channels and the number of output channels
% can be different, for example when the user specifies that he/she wants
% to add the implicit EEG reference channel to the data matrix.
%
% The filtering of the data can introduce artifacts at the edges, hence it
% is better to pad the data with some extra signal at the begin and end.
% After filtering, this padding is removed and the other preprocessing
% steps are applied to the remainder of the data. The input fields
% begpadding and endpadding should be specified in samples. You can also
% leave them empty, which implies that the data is not padded.
%
% The configuration can contain
%   cfg.lpfilter      = 'no' or 'yes'  lowpass filter
%   cfg.hpfilter      = 'no' or 'yes'  highpass filter
%   cfg.bpfilter      = 'no' or 'yes'  bandpass filter
%   cfg.bsfilter      = 'no' or 'yes'  bandstop filter
%   cfg.dftfilter     = 'no' or 'yes'  line noise removal using discrete fourier transform
%   cfg.medianfilter  = 'no' or 'yes'  jump preserving median filter
%   cfg.lpfreq        = lowpass  frequency in Hz
%   cfg.hpfreq        = highpass frequency in Hz
%   cfg.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%   cfg.bsfreq        = bandstop frequency range, specified as [low high] in Hz
%   cfg.dftfreq       = line noise frequencies for DFT filter, default [50 100 150] Hz
%   cfg.lpfiltord     = lowpass  filter order (default set in low-level function)
%   cfg.hpfiltord     = highpass filter order (default set in low-level function)
%   cfg.bpfiltord     = bandpass filter order (default set in low-level function)
%   cfg.bsfiltord     = bandstop filter order (default set in low-level function)
%   cfg.medianfiltord = length of median filter
%   cfg.lpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.hpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.bpfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.bsfilttype    = digital filter type, 'but' (default) or 'firws' or 'fir' or 'firls'
%   cfg.lpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.hpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.bpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.bsfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
%   cfg.lpinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.hpinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.bpinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.bsinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
%   cfg.lpfiltdf      = lowpass transition width (firws, overrides order, default set in low-level function)
%   cfg.hpfiltdf      = highpass transition width (firws, overrides order, default set in low-level function)
%   cfg.bpfiltdf      = bandpass transition width (firws, overrides order, default set in low-level function)
%   cfg.bsfiltdf      = bandstop transition width (firws, overrides order, default set in low-level function)
%   cfg.lpfiltwintype = lowpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.hpfiltwintype = highpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.bpfiltwintype = bandpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.bsfiltwintype = bandstop window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
%   cfg.lpfiltdev     = lowpass max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.hpfiltdev     = highpass max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.bpfiltdev     = bandpass max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.bsfiltdev     = bandstop max passband deviation (firws with 'kaiser' window, default 0.001 set in low-level function)
%   cfg.dftreplace    = 'zero' or 'neighbour', method used to reduce line noise, 'zero' implies DFT filter, 'neighbour' implies spectrum interpolation (default = 'zero')
%   cfg.dftbandwidth  = bandwidth of line noise frequencies, applies to spectrum interpolation, in Hz (default = [1 2 3])
%   cfg.dftneighbourwidth = bandwidth of frequencies neighbouring line noise frequencies, applies to spectrum interpolation, in Hz (default = [2 2 2])
%   cfg.plotfiltresp  = 'no' or 'yes', plot filter responses (firws, default = 'no')
%   cfg.usefftfilt    = 'no' or 'yes', use fftfilt instead of filter (firws, default = 'no')
%   cfg.demean        = 'no' or 'yes'
%   cfg.baselinewindow = [begin end] in seconds, the default is the complete trial
%   cfg.detrend       = 'no' or 'yes', this is done on the complete trial
%   cfg.polyremoval   = 'no' or 'yes', this is done on the complete trial
%   cfg.polyorder     = polynome order (default = 2)
%   cfg.derivative    = 'no' (default) or 'yes', computes the first order derivative of the data
%   cfg.hilbert       = 'no', 'abs', 'complex', 'real', 'imag', 'absreal', 'absimag' or 'angle' (default = 'no')
%   cfg.rectify       = 'no' or 'yes'
%   cfg.precision     = 'single' or 'double' (default = 'double')
%   cfg.absdiff       = 'no' or 'yes', computes absolute derivative (i.e.first derivative then rectify)
%
% Preprocessing options that you should only use for EEG data are
%   cfg.reref         = 'no' or 'yes' (default = 'no')
%   cfg.refchannel    = cell-array with new EEG reference channel(s)
%   cfg.refmethod     = 'avg', 'median', or 'bipolar' (default = 'avg')
%   cfg.implicitref   = 'label' or empty, add the implicit EEG reference as zeros (default = [])
%   cfg.montage       = 'no' or a montage structure (default = 'no')
%
% See also FT_READ_DATA, FT_READ_HEADER

% TODO implement decimation and/or resampling

% Copyright (C) 2004-2012, Robert Oostenveld
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

% compute fsample
fsample = 1./nanmean(diff(ttime));

if nargin<5 || isempty(begpadding)
    begpadding = 0;
end
if nargin<6 || isempty(endpadding)
    endpadding = 0;
end

if iscell(cfg)
    % recurse over the subsequent preprocessing stages
    if begpadding>0 || endpadding>0
        ft_error('multiple preprocessing stages are not supported in combination with filter padding');
    end
    for i=1:length(cfg)
        tmpcfg = cfg{i};
        if nargout==1
            [dat                     ] = preproc(dat, label, ttime, tmpcfg, begpadding, endpadding);
        elseif nargout==2
            [dat, label              ] = preproc(dat, label, ttime, tmpcfg, begpadding, endpadding);
        elseif nargout==3
            [dat, label, ttime        ] = preproc(dat, label, ttime, tmpcfg, begpadding, endpadding);
        elseif nargout==4
            [dat, label, ttime, tmpcfg] = preproc(dat, label, ttime, tmpcfg, begpadding, endpadding);
            cfg{i} = tmpcfg;
        end
    end
    % ready with recursing over the subsequent preprocessing stages
    return
end

% set the defaults for the rereferencing options
cfg.reref =                ft_getopt(cfg, 'reref', 'no');
cfg.refchannel =           ft_getopt(cfg, 'refchannel', {});
cfg.refmethod =            ft_getopt(cfg, 'refmethod', 'avg');
cfg.implicitref =          ft_getopt(cfg, 'implicitref', []);
% set the defaults for the signal processing options
cfg.polyremoval =          ft_getopt(cfg, 'polyremoval', 'no');
cfg.polyorder =            ft_getopt(cfg, 'polyorder', 2);
cfg.detrend =              ft_getopt(cfg, 'detrend', 'no');
cfg.demean =               ft_getopt(cfg, 'demean', 'no');
cfg.baselinewindow =       ft_getopt(cfg, 'baselinewindow', 'all');
cfg.dftfilter =            ft_getopt(cfg, 'dftfilter', 'no');
cfg.lpfilter =             ft_getopt(cfg, 'lpfilter', 'no');
cfg.hpfilter =             ft_getopt(cfg, 'hpfilter', 'no');
cfg.bpfilter =             ft_getopt(cfg, 'bpfilter', 'no');
cfg.bsfilter =             ft_getopt(cfg, 'bsfilter', 'no');
cfg.lpfiltord =            ft_getopt(cfg, 'lpfiltord', []);
cfg.hpfiltord =            ft_getopt(cfg, 'hpfiltord', []);
cfg.bpfiltord =            ft_getopt(cfg, 'bpfiltord', []);
cfg.bsfiltord =            ft_getopt(cfg, 'bsfiltord', []);
cfg.lpfilttype =           ft_getopt(cfg, 'lpfilttype', 'but');
cfg.hpfilttype =           ft_getopt(cfg, 'hpfilttype', 'but');
cfg.bpfilttype =           ft_getopt(cfg, 'bpfilttype', 'but');
cfg.bsfilttype =           ft_getopt(cfg, 'bsfilttype', 'but');
if strcmp(cfg.lpfilttype, 'firws'), cfg.lpfiltdir = ft_getopt(cfg, 'lpfiltdir', 'onepass-zerophase'); else, cfg.lpfiltdir = ft_getopt(cfg, 'lpfiltdir', 'twopass'); end
if strcmp(cfg.hpfilttype, 'firws'), cfg.hpfiltdir = ft_getopt(cfg, 'hpfiltdir', 'onepass-zerophase'); else, cfg.hpfiltdir = ft_getopt(cfg, 'hpfiltdir', 'twopass'); end
if strcmp(cfg.bpfilttype, 'firws'), cfg.bpfiltdir = ft_getopt(cfg, 'bpfiltdir', 'onepass-zerophase'); else, cfg.bpfiltdir = ft_getopt(cfg, 'bpfiltdir', 'twopass'); end
if strcmp(cfg.bsfilttype, 'firws'), cfg.bsfiltdir = ft_getopt(cfg, 'bsfiltdir', 'onepass-zerophase'); else, cfg.bsfiltdir = ft_getopt(cfg, 'bsfiltdir', 'twopass'); end
cfg.lpinstabilityfix =     ft_getopt(cfg, 'lpinstabilityfix', 'no');
cfg.hpinstabilityfix =     ft_getopt(cfg, 'hpinstabilityfix', 'no');
cfg.bpinstabilityfix =     ft_getopt(cfg, 'bpinstabilityfix', 'no');
cfg.bsinstabilityfix =     ft_getopt(cfg, 'bsinstabilityfix', 'no');
cfg.lpfiltdf =             ft_getopt(cfg, 'lpfiltdf',  []);
cfg.hpfiltdf =             ft_getopt(cfg, 'hpfiltdf', []);
cfg.bpfiltdf =             ft_getopt(cfg, 'bpfiltdf', []);
cfg.bsfiltdf =             ft_getopt(cfg, 'bsfiltdf', []);
cfg.lpfiltwintype =        ft_getopt(cfg, 'lpfiltwintype', 'hamming');
cfg.hpfiltwintype =        ft_getopt(cfg, 'hpfiltwintype', 'hamming');
cfg.bpfiltwintype =        ft_getopt(cfg, 'bpfiltwintype', 'hamming');
cfg.bsfiltwintype =        ft_getopt(cfg, 'bsfiltwintype', 'hamming');
cfg.lpfiltdev =            ft_getopt(cfg, 'lpfiltdev', []);
cfg.hpfiltdev =            ft_getopt(cfg, 'hpfiltdev', []);
cfg.bpfiltdev =            ft_getopt(cfg, 'bpfiltdev', []);
cfg.bsfiltdev =            ft_getopt(cfg, 'bsfiltdev', []);
cfg.plotfiltresp =         ft_getopt(cfg, 'plotfiltresp', 'no');
cfg.usefftfilt =           ft_getopt(cfg, 'usefftfilt', 'no');
cfg.medianfilter =         ft_getopt(cfg, 'medianfilter ', 'no');
cfg.medianfiltord =        ft_getopt(cfg, 'medianfiltord', 9);
cfg.dftfreq =              ft_getopt(cfg, 'dftfreq', [50 100 150]);
cfg.hilbert =              ft_getopt(cfg, 'hilbert', 'no');
cfg.derivative =           ft_getopt(cfg, 'derivative', 'no');
cfg.rectify =              ft_getopt(cfg, 'rectify', 'no');
cfg.boxcar =               ft_getopt(cfg, 'boxcar', 'no');
cfg.absdiff =              ft_getopt(cfg, 'absdiff', 'no');
cfg.precision =            ft_getopt(cfg, 'precision', []);
cfg.conv =                 ft_getopt(cfg, 'conv', 'no');
cfg.montage =              ft_getopt(cfg, 'montage', 'no');
cfg.dftinvert =            ft_getopt(cfg, 'dftinvert', 'no');
cfg.standardize =          ft_getopt(cfg, 'standardize', 'no');
cfg.denoise =              ft_getopt(cfg, 'denoise', '');
cfg.subspace =             ft_getopt(cfg, 'subspace', []);
cfg.custom =               ft_getopt(cfg, 'custom', '');
cfg.resample =             ft_getopt(cfg, 'resample', '');

% test whether the MATLAB signal processing toolbox is available
if strcmp(cfg.medianfilter, 'yes') && ~ft_hastoolbox('signal')
    ft_error('median filtering requires the MATLAB signal processing toolbox');
end

% do a sanity check on the filter configuration
if strcmp(cfg.bpfilter, 'yes') && ...
        (strcmp(cfg.hpfilter, 'yes') || strcmp(cfg.lpfilter,'yes'))
    ft_error('you should not apply both a bandpass AND a lowpass/highpass filter');
end

% do a sanity check on the hilbert transform configuration
if strcmp(cfg.hilbert, 'yes') && ~strcmp(cfg.bpfilter, 'yes')
    ft_warning('Hilbert transform should be applied in conjunction with bandpass filter')
end

% do a sanity check on hilbert and rectification
if strcmp(cfg.hilbert, 'yes') && strcmp(cfg.rectify, 'yes')
    ft_error('Hilbert transform and rectification should not be applied both')
end

% do a sanity check on the rereferencing/montage
if ~strcmp(cfg.reref, 'no') && ~strcmp(cfg.montage, 'no')
    ft_error('cfg.reref and cfg.montage are mutually exclusive')
end

% lnfilter is no longer used
if isfield(cfg, 'lnfilter') && strcmp(cfg.lnfilter, 'yes')
    ft_error('line noise filtering using the option cfg.lnfilter is not supported any more, use cfg.bsfilter instead')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the rereferencing in case of EEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(cfg.implicitref) && ~any(match_str(cfg.implicitref,label))
    label = [label(:)' {cfg.implicitref}]';
    dat(end+1,:) = 0;
end

if strcmp(cfg.reref, 'yes')
    if strcmp(cfg.refmethod, 'bipolar')
        % this is implemented as a montage that the user does not get to see
        tmpcfg = keepfields(cfg, {'refmethod', 'implicitref', 'refchannel', 'channel'});
        tmpcfg.showcallinfo = 'no';
        montage = ft_prepare_montage(tmpcfg);
        % convert the data temporarily to a raw structure
        tmpdata.trial = {dat};
        tmpdata.time  = {time};
        tmpdata.label = label;
        % apply the montage to the data
        tmpdata = ft_apply_montage(tmpdata, montage, 'feedback', 'none');
        dat   = tmpdata.trial{1}; % the number of channels can have changed
        label = tmpdata.label;    % the output channels can be different than the input channels
        clear tmpdata
    else
        % mean or median based derivation of specified or all channels
        cfg.refchannel = ft_channelselection(cfg.refchannel, label);
        refindx = match_str(label, cfg.refchannel);
        if isempty(refindx)
            ft_error('reference channel was not found')
        end
        dat = ft_preproc_rereference(dat, refindx, cfg.refmethod);
    end
end

if ~strcmp(cfg.montage, 'no') && ~isempty(cfg.montage)
    % convert the data temporarily to a raw structure
    tmpdata.trial = {dat};
    tmpdata.time  = {ttime};
    tmpdata.label = label;
    % apply the montage to the data
    tmpdata = ft_apply_montage(tmpdata, cfg.montage, 'feedback', 'none');
    dat   = tmpdata.trial{1}; % the number of channels can have changed
    label = tmpdata.label;    % the output channels can be different than the input channels
    clear tmpdata
end

if any(any(isnan(dat)))
    % filtering is not possible for at least a selection of the data
    ft_warning('data contains NaNs, no filtering or preprocessing applied');
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % do the filtering on the padded data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(cfg.denoise)
        hflag    = isfield(cfg.denoise, 'hilbert') && strcmp(cfg.denoise.hilbert, 'yes');
        datlabel = match_str(label, cfg.denoise.channel);
        reflabel = match_str(label, cfg.denoise.refchannel);
        tmpdat   = ft_preproc_denoise(dat(datlabel,:), dat(reflabel,:), hflag);
        dat(datlabel,:) = tmpdat;
    end
    
    % The filtering should in principle be done prior to the demeaning to
    % ensure that the resulting mean over the baseline window will be
    % guaranteed to be zero (even if there are filter artifacts).
    % However, the filtering benefits from the data being pulled towards zero,
    % causing less edge artifacts. That is why we start by removing the slow
    % drift, then filter, and then repeat the demean/detrend/polyremove.
    if strcmp(cfg.polyremoval, 'yes')
        nsamples  = size(dat,2);
        begsample = 1        + begpadding;
        endsample = nsamples - endpadding;
        dat = ft_preproc_polyremoval(dat, cfg.polyorder, begsample, endsample); % this will also demean and detrend
    elseif strcmp(cfg.detrend, 'yes')
        nsamples  = size(dat,2);
        begsample = 1        + begpadding;
        endsample = nsamples - endpadding;
        dat = ft_preproc_polyremoval(dat, 1, begsample, endsample); % this will also demean
    elseif strcmp(cfg.demean, 'yes')
        nsamples  = size(dat,2);
        begsample = 1        + begpadding;
        endsample = nsamples - endpadding;
        dat = ft_preproc_polyremoval(dat, 0, begsample, endsample);
    end
    
    if strcmp(cfg.medianfilter, 'yes'), dat = ft_preproc_medianfilter(dat, cfg.medianfiltord); end
    if strcmp(cfg.lpfilter, 'yes'),     dat = ft_preproc_lowpassfilter(dat, fsample, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype, cfg.lpfiltdir, cfg.lpinstabilityfix, cfg.lpfiltdf, cfg.lpfiltwintype, cfg.lpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
    if strcmp(cfg.hpfilter, 'yes'),     dat = ft_preproc_highpassfilter(dat, fsample, cfg.hpfreq, cfg.hpfiltord, cfg.hpfilttype, cfg.hpfiltdir, cfg.hpinstabilityfix, cfg.hpfiltdf, cfg.hpfiltwintype, cfg.hpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
    if strcmp(cfg.bpfilter, 'yes'),     dat = ft_preproc_bandpassfilter(dat, fsample, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype, cfg.bpfiltdir, cfg.bpinstabilityfix, cfg.bpfiltdf, cfg.bpfiltwintype, cfg.bpfiltdev, cfg.plotfiltresp, cfg.usefftfilt); end
    if strcmp(cfg.bsfilter, 'yes')
        for i=1:size(cfg.bsfreq,1)
            % apply a bandstop filter for each of the specified bands, i.e. cfg.bsfreq should be Nx2
            dat = ft_preproc_bandstopfilter(dat, fsample, cfg.bsfreq(i,:), cfg.bsfiltord, cfg.bsfilttype, cfg.bsfiltdir, cfg.bsinstabilityfix, cfg.bsfiltdf, cfg.bsfiltwintype, cfg.bsfiltdev, cfg.plotfiltresp, cfg.usefftfilt);
        end
    end
    if strcmp(cfg.polyremoval, 'yes')
        % the begin and endsample of the polyremoval period correspond to the complete data minus padding
        nsamples  = size(dat,2);
        begsample = 1        + begpadding;
        endsample = nsamples - endpadding;
        dat = ft_preproc_polyremoval(dat, cfg.polyorder, begsample, endsample);
    end
    if strcmp(cfg.detrend, 'yes')
        % the begin and endsample of the detrend period correspond to the complete data minus padding
        nsamples  = size(dat,2);
        begsample = 1        + begpadding;
        endsample = nsamples - endpadding;
        dat = ft_preproc_detrend(dat, begsample, endsample);
    end
    if strcmp(cfg.demean, 'yes')
        if ischar(cfg.baselinewindow) && strcmp(cfg.baselinewindow, 'all')
            % the begin and endsample of the baseline period correspond to the complete data minus padding
            nsamples  = size(dat,2);
            begsample = 1        + begpadding;
            endsample = nsamples - endpadding;
            dat       = ft_preproc_baselinecorrect(dat, begsample, endsample);
        else
            % determine the begin and endsample of the baseline period and baseline correct for it
            begsample = nearest(ttime, cfg.baselinewindow(1));
            endsample = nearest(ttime, cfg.baselinewindow(2));
            dat       = ft_preproc_baselinecorrect(dat, begsample, endsample);
        end
    end
    if strcmp(cfg.dftfilter, 'yes')
        datorig = dat;
        optarg = {};
        if isfield(cfg, 'dftreplace')
            optarg = cat(2, optarg, {'dftreplace', cfg.dftreplace});
            if strcmp(cfg.dftreplace, 'neighbour') && (begpadding>0 || endpadding>0)
                ft_error('Padding by data mirroring is not supported for spectrum interpolation.');
            end
        end
        if isfield(cfg, 'dftbandwidth')
            optarg = cat(2, optarg, {'dftbandwidth', cfg.dftbandwidth});
        end
        if isfield(cfg, 'dftneighbourwidth')
            optarg = cat(2, optarg, {'dftneighbourwidth', cfg.dftneighbourwidth});
        end
        dat     = ft_preproc_dftfilter(dat, fsample, cfg.dftfreq, optarg{:});
        if strcmp(cfg.dftinvert, 'yes')
            dat = datorig - dat;
        end
    end
    if ~strcmp(cfg.hilbert, 'no')
        dat = ft_preproc_hilbert(dat, cfg.hilbert);
    end
    if strcmp(cfg.rectify, 'yes')
        dat = ft_preproc_rectify(dat);
    end
    if isnumeric(cfg.boxcar)
        numsmp = round(cfg.boxcar*fsample);
        if ~rem(numsmp,2)
            % the kernel should have an odd number of samples
            numsmp = numsmp+1;
        end
        % kernel = ones(1,numsmp) ./ numsmp;
        % dat    = convn(dat, kernel, 'same');
        dat = ft_preproc_smooth(dat, numsmp); % better edge behavior
    end
    if isnumeric(cfg.conv)
        kernel = (cfg.conv(:)'./sum(cfg.conv));
        if ~rem(length(kernel),2)
            kernel = [kernel 0];
        end
        dat = convn(dat, kernel, 'same');
    end
    if strcmp(cfg.derivative, 'yes')
        dat = ft_preproc_derivative(dat, 1);
    end
    if strcmp(cfg.absdiff, 'yes')
        % this implements abs(diff(data), which is required for jump detection
        dat = abs([diff(dat, 1, 2) zeros(size(dat,1),1)]);
    end
    if strcmp(cfg.standardize, 'yes')
        dat = ft_preproc_standardize(dat, 1, size(dat,2));
    end
    if ~isempty(cfg.subspace)
        dat = ft_preproc_subspace(dat, cfg.subspace);
    end
    if ~isempty(cfg.custom)
        if ~isfield(cfg.custom, 'nargout')
            cfg.custom.nargout = 1;
        end
        if cfg.custom.nargout==1
            dat = feval(cfg.custom.funhandle, dat, cfg.custom.varargin);
        elseif cfg.custom.nargout==2
            [dat, ttime] = feval(cfg.custom.funhandle, dat, cfg.custom.varargin);
        end
    end
    if strcmp(cfg.resample, 'yes')
        if ~isfield(cfg, 'resamplefs')
            cfg.resamplefs = fsample./2;
        end
        if ~isfield(cfg, 'resamplemethod')
            cfg.resamplemethod = 'resample';
        end
        [dat               ] = ft_preproc_resample(dat,  fsample, cfg.resamplefs, cfg.resamplemethod);
        [ttime, ~, ~] = ft_preproc_resample(ttime, fsample, cfg.resamplefs, cfg.resamplemethod);
    end
    if ~isempty(cfg.precision)
        % convert the data to another numeric precision, i.e. double, single or int32
        dat = cast(dat, cfg.precision);
    end
end % if any(isnan)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the filter padding and do the preprocessing on the remaining trial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if begpadding~=0 || endpadding~=0
    dat = ft_preproc_padding(dat, 'remove', begpadding, endpadding);
    if strcmp(cfg.demean, 'yes') || nargout>2
        ttime = ft_preproc_padding(ttime, 'remove', begpadding, endpadding);
    end
end
end

function [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = do_plot_chantrl(info,level)


[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = set_maxper(level, info.chansel, info.trlsel, strcmp(info.metric, 'min'));

if info.pflag == 1
    figure,
    axis ij;
    ymax = max(maxperchan); ymin = min(maxperchan); xmax = info.nchan;
    subplot 121,
    plot(maxperchan,'.'), xlabel('Channel number'), ylabel(info.metric),
    axis([0.5 xmax+0.5 0.8*ymin 1.2*ymax]);
    title('Channel')
    subplot 122,
    plot(maxpertrl,'.'), xlabel('Trial number'), ylabel(info.metric);
    xmax = info.ntrl; ymax = max(maxpertrl); ymin = min(maxpertrl);
    axis([0.5 xmax+0.5 (1-sign(ymin)*0.2)*ymin (1+sign(ymax)*0.2)*ymax]);
    title('Trial')
    sgtitle(info.metric)
end


    function [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = set_maxper(level, chansel, trlsel, minflag)
        if minflag
            % take the negative maximum, i.e. the minimum
            level = -1 * level;
        end
        % determine the maximum value
        maxperchan_all = max(level, [], 2);
        maxpertrl_all  = max(level, [], 1);
        % determine the maximum value over the remaining selection
        level(~chansel, :) = nan;
        level(:, ~trlsel)  = nan;
        maxperchan     = max(level, [], 2);
        maxpertrl      = max(level, [], 1);
        if minflag
            maxperchan     = -1 * maxperchan;
            maxpertrl      = -1 * maxpertrl;
            maxperchan_all = -1 * maxperchan_all;
            maxpertrl_all  = -1 * maxpertrl_all;
            level          = -1 * level;
        end
    end
end
