function varargout = process_ft_ica_clean(varargin )
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
sProcess.Comment     = 'FieldTrip: process_ft_ica_clean';
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

cfg            = [];
cfg.method     = 'runica';
cfg.numcomponent = 20;       % specify the component(s) that should be plotted
comp           = ft_componentanalysis(cfg, ftData1);

%%
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);

cfg = [];
cfg.viewmode = 'component';
cfg.layout = lay;
ft_databrowser(cfg, comp);

%%
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';%compute the power spectrum in all ICs
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = 2:2:30;
freq = ft_freqanalysis(cfg, comp);

%%
n = 20;
nby1 = 5; nby2 = 4;

Nfigs = ceil(size(comp.topo,1)/n);
tot = Nfigs*n;

rptvect = 1:size(comp.topo,1);
rptvect = padarray(rptvect, [0 tot-size(comp.topo,1)], 0,'post');
rptvect = reshape(rptvect,n,Nfigs)';

figure
for r=1:n
    cfg = [];
    cfg.channel = rptvect(:,r);
    subplot(nby1,nby2,r);
    set(gca,'color','none');
    plot(freq.freq, freq.powspctrm(r,:))
    xlim([min(freq.freq), max(freq.freq)])
    title(['IC', num2str(r)])
end

%%
disp('Enter bad ICs (use [] for no correction):')
bic = input(['']);

if ~ isempty(bic)
    cfg = [];
    cfg.component = comp.label(bic);
    cfg.updatesens = 'no';
    ic_data = ft_rejectcomponent(cfg, comp, ftData1);
    
    %%
    ProtocolInfo = bst_get('ProtocolInfo');
    [datapath,~] = fileparts(DataFile);
    [~,srcPath] = fileparts(datapath);
    srcDir  = bst_fullfile(ProtocolInfo.STUDIES, [srcPath,'_ica']);
    
    cd(fullfile(ProtocolInfo.STUDIES, datapath))
    for iInput = 1:length(sInputs)
        DataFile = sInputs(iInput).FileName;
        D = load(file_fullpath(DataFile));
        [a,FileName] = fileparts(DataFile);
        FileName_new = [FileName, '_ic'];
        tkz = tokenize(D.Comment, ' ');
        D.Comment = [tkz{1} '_ica ',tkz{2}];
        D.F(iChannelsData,:) = ic_data.trial{iInput};
        save(FileName_new, '-struct', 'D');
        
    end
    disp('done, reload datafile!')
    
else,
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
