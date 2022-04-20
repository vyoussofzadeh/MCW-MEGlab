function varargout = process_ft_sourceanalysis_MNE_ICA(varargin )
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
sProcess.Comment     = 'FieldTrip: ft_sourceanalysis MNE ICA v052521';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 357;
sProcess.Description = 'https://github.com/vyoussofzadeh/MCW-MEGlab';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data'};
sProcess.OutputTypes = {'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

% Label: Time
sProcess.options.label1.Comment = '<BR><B>Time of interest:</B>';
sProcess.options.label1.Type    = 'label';
% Active time window
sProcess.options.poststim.Comment = 'Time interval:';
sProcess.options.poststim.Type    = 'poststim';
sProcess.options.poststim.Value   = [];

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
% Inverse options
PostStim = sProcess.options.poststim.Value{1};
Modality = sProcess.options.sensortype.Value{1};
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

% ===== FIELDTRIP: ft_freqanalysis =====
bst_progress('text', 'Calling FieldTrip function: ft_freqanalysis...');

%%
TT = ftData.time;

for i=1:length(TT)
    TT{i} = ftData.time{1};
end
ftData.time = TT;

%%
cfg = [];
cfg.toilim = PostStim;
datain = ft_redefinetrial(cfg, ftData);


disp('hpfreq to lpfreq eg, [18,23] in Hz:');
disp('[] for no filtering:');
foi = input('');

if ~isempty(foi)
    cfg        = [];
    cfg.demean = 'yes';
    cfg.dftfilter = 'yes';
    cfg.hpfilter = 'yes';
    cfg.lpfilter = 'yes';
    cfg.hpfiltord = 3;
    cfg.lpfreq = foi(2);
    cfg.hpfreq = foi(1);
    %         cfg.channel = {'megmag', 'meggrad'};
    fcln_data        = ft_preprocessing(cfg, datain);
else
    fcln_data = datain;
end

%%
disp('number ICs to be extracted (e.g., 20):');
nic = input('');

%%
ftData1 = [];
ftData1.label = fcln_data.label;
%         ftData1.grad = ftData.grad;
ftData1.dimord = 'chan_time';
ftData1.trial = fcln_data.trial;
ftData1.time = fcln_data.time;

%%
cfg            = [];
cfg.method     = 'runica';
cfg.numcomponent = nic;       % specify the component(s) that should be plotted
% cfg.numcomponent = 1;       % specify the component(s) that should be plotted
comp           = ft_componentanalysis(cfg, ftData1);

%%
OutputDir = bst_fileparts(file_fullpath(DataFile));
Index = strfind(OutputDir, 'data/');
bsdir = OutputDir(1:Index(end)-1);
bsdatadir = fullfile(bsdir,'data');
bsanatdir = fullfile(bsdir,'anat');


%% Step8: head (forward) model
sourcemodel = ft_read_headshape(fullfile(bsanatdir,HeadModelMat.SurfaceFile));
[ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);

%%
disp('Fast post-pre contrast: 1, slow cluster-based statistics: 2 ?');
%         st = input('');
st = 1;
if st == 2, Method = 'dics_stat'; end

Method = 'mne';

%% step9: Surface-based source analysis
switch Method
    case 'mne'
        
        %%
        % Prepare the component data in order for ft_sourceanalysis to be able to
        % swallow it
        mixing   = comp.topo;
        channels = comp.topolabel;
        % normalisation of the topographies
        for i = 1:size(mixing, 2)
            val(i) = 0.01*max(abs(mixing(:, i)));
            mixing(:, i) = mixing(:, i)/val(i);
        end
        
        % create a 'timelock' structure
        tlck = [];
        tlck.label = channels;
        tlck.cov = eye(numel(tlck.label));
        tlck.time=1;
        tlck.grad = ftData.grad;
        tlck.dimord = 'chan_time';
        
        
        %%
        % do an MNE with different regularisation for each component
        
        % this parameter is hard-coded here
        noise_level = 8;
        
        % specify the static part of the cfg for the source reconstruction
        cfg               = [];
        cfg.method        = 'mne';
        % cfg.grid          = gridLF;
        % cfg.vol           = individual_headmodel;
        
        cfg.sourcemodel = ftLeadfield;
        cfg.headmodel = ftHeadmodel;
        
        %         cfg.grid = individual_grid;
        %         cfg.headmodel = individual_headmodel;
        cfg.channel       = channels;
        cfg.mne.prewhiten = 'yes';
        cfg.mne.noisecov  = eye(numel(channels))*noise_level;
        
        % loop over components, due to component-specific regularisation
        for i=1:size(mixing,2)
            
            % use the channel-level topography of the current component
            tlck.avg = mixing(:,i);
            
            % estimate the snr of the current component
            cfg.mne.snr = sqrt(mean((mixing(:,i)-mean(mixing(:,i))).^2))/noise_level;
            noisevec(i) = cfg.mne.snr;
            
            tmp = ft_sourceanalysis(cfg, tlck);
            inside_indices = find(tmp.inside(:))';
            if i==1
                % create the output source structure in the first iteration
                source=tmp;
            else
                % concatenate the reconstructed source level topography to the previously computed ones
                for k = 1:numel(inside_indices)
                    source.avg.mom{inside_indices(k)} = cat(2,source.avg.mom{inside_indices(k)}, tmp.avg.mom{inside_indices(k)});
                end
                source.avg.pow = horzcat(source.avg.pow,tmp.avg.pow);
            end
        end
        
        %%
        cfg = [];
        cfg.layout = 'neuromag306mag.lay';
        lay = ft_prepare_layout(cfg);
        
        
        % pause,
        close all,
        cfg = [];
        cfg.viewmode = 'component';
        cfg.layout = lay;
        ft_databrowser(cfg, comp);
        %                 set(gcf, 'Position', [200   400   1200   800]);
        %                 colormap(brewermap(256, '*RdYlBu'));
        %         title(subj)
        
        %%
        switch HeadModelMat.HeadModelType
            case 'surface'
                disp('enter ICs:')
                nIC = input(' ');
                for i=1:length(nIC)
                    source1 = [];
                    source1.pow = source.avg.pow(:,nIC(i));
                    
                    figure
                    m = source1.pow;
                    bnd.pnt = sourcemodel.pos;
                    bnd.tri = sourcemodel.tri;
                    ft_plot_mesh(bnd, 'vertexcolor', abs(m));
                    colorbar
                    view([-180,0])
                    title(['IC:', num2str(nIC(i))])
                end
        end
        
end
disp('enter ICs for saving surfaces:')
nIC = input(' ');

disp('avg only, yes(avg of selected ICs only)=1, no=2')
avg_ask = input(' ');

close all,

ResultsMat = db_template('resultsmat');

%%
bst_progress('text', 'Saving source file...');
bst_progress('inc', 1);
% Output study
if (length(sInputs) == 1)
    iStudyOut = sInputs(1).iStudy;
    RefDataFile = sInputs(iChanInputs(iInput)).FileName;
else
    [tmp, iStudyOut] = bst_process('GetOutputStudy', sProcess, sInputs);
    RefDataFile = [];
end
%%

if avg_ask == 1
    
    ResultsMat.ImagingKernel = [];
    switch Method
        case 'mne'
            ResultsMat.ImageGridAmp  = mean(source.avg.pow(:,nIC),2);
            ResultsMat.cfg           = source.cfg;
    end
    ResultsMat.nComponents   = 1;
    ResultsMat.Function      = Method;
    ResultsMat.Time          = 1;
    ResultsMat.DataFile      = DataFile;
    ResultsMat.HeadModelFile = HeadModelFile;
    ResultsMat.HeadModelType = HeadModelMat.HeadModelType;
    ResultsMat.ChannelFlag   = DataMat.ChannelFlag;
    ResultsMat.GoodChannel   = iChannelsData;
    ResultsMat.SurfaceFile   = HeadModelMat.SurfaceFile;
    ResultsMat.nAvg          = DataMat.nAvg;
    ResultsMat.Leff          = DataMat.Leff;
    ResultsMat.Comment       = ['ft_sourceanalysis:' Method, '_IC_avg_', num2str(length(nIC)), '_', num2str(foi(1)),'_',num2str(foi(2)),'_Hz_',datestr(now, 'dd/mm/yy-HH:MM')];
    switch lower(ResultsMat.HeadModelType)
        case 'volume'
            ResultsMat.GridLoc    = HeadModelMat.GridLoc;
            % ResultsMat.GridOrient = [];
        case 'surface'
            ResultsMat.GridLoc    = [];
            % ResultsMat.GridOrient = [];
        case 'mixed'
            ResultsMat.GridLoc    = HeadModelMat.GridLoc;
            ResultsMat.GridOrient = HeadModelMat.GridOrient;
    end
    ResultsMat = bst_history('add', ResultsMat, 'compute', ['ft_sourceanalysis: ' Method ' ' Modality]);
    
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
else
    
    for i=1:length(nIC)
        
        source1 = [];
        source1.pow = source.avg.pow(:,nIC(i));
        
        % === CREATE OUTPUT STRUCTURE ===
        bst_progress('text', 'Saving source file...');
        bst_progress('inc', 1);
        % Create structure
        %             ResultsMat = db_template('resultsmat');
        ResultsMat.ImagingKernel = [];
        switch Method
            case 'lcmv'
                ResultsMat.ImageGridAmp  = abs((source.(mask)));
                ResultsMat.cfg           = source.cfg;
            case 'dics_stat'
                ResultsMat.ImageGridAmp  = abs((stats2.stat));
                ResultsMat.cfg           = stat.cfg;
            case 'mne'
                ResultsMat.ImageGridAmp  = source1.pow;
                ResultsMat.cfg           = source.cfg;
        end
        ResultsMat.nComponents   = 1;
        ResultsMat.Function      = Method;
        ResultsMat.Time          = 1;
        ResultsMat.DataFile      = DataFile;
        ResultsMat.HeadModelFile = HeadModelFile;
        ResultsMat.HeadModelType = HeadModelMat.HeadModelType;
        ResultsMat.ChannelFlag   = DataMat.ChannelFlag;
        ResultsMat.GoodChannel   = iChannelsData;
        ResultsMat.SurfaceFile   = HeadModelMat.SurfaceFile;
        ResultsMat.nAvg          = DataMat.nAvg;
        ResultsMat.Leff          = DataMat.Leff;
        ResultsMat.Comment       = ['ft_sourceanalysis:' Method, '_IC', num2str(nIC(i)), '_', num2str(foi(1)),'_',num2str(foi(2)),'_Hz_',datestr(now, 'dd/mm/yy-HH:MM')];
        %         ResultsMat.Comment       = ['ft_sourceanalysis: ' Method, '_',num2str(foi(1)),'_',num2str(foi(2)),'_Hz_', num2str(toi(2,1)),'_',num2str(toi(2,2)), 'sec ', datestr(now, 'dd/mm/yy-HH:MM')];
        switch lower(ResultsMat.HeadModelType)
            case 'volume'
                ResultsMat.GridLoc    = HeadModelMat.GridLoc;
                % ResultsMat.GridOrient = [];
            case 'surface'
                ResultsMat.GridLoc    = [];
                % ResultsMat.GridOrient = [];
            case 'mixed'
                ResultsMat.GridLoc    = HeadModelMat.GridLoc;
                ResultsMat.GridOrient = HeadModelMat.GridOrient;
        end
        ResultsMat = bst_history('add', ResultsMat, 'compute', ['ft_sourceanalysis: ' Method ' ' Modality]);
        
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
    end
end
% Save database
db_save();
% Hide progress bar
bst_progress('stop');
end


function [time_of_interest,freq_of_interest] = vy_tfr_plot(cfg_main, tfr)

%%
% First compute the average over trials:
cfg = [];
freq_avg = ft_freqdescriptives(cfg, tfr);

if cfg_main.bslcorr == 1
    % And baseline-correct the average:
    cfg = [];
    cfg.baseline = [-0.3 0];
    cfg.baselinetype = 'db'; % Use decibel contrast here
    freq_avg_bsl = ft_freqbaseline(cfg, freq_avg);
    freq_avg_bsl.powspctrm(isnan(freq_avg_bsl.powspctrm))=0;
    meanpow = squeeze(mean(freq_avg_bsl.powspctrm, 1));
    
else
    freq_avg.powspctrm(isnan(freq_avg.powspctrm))=0;
    meanpow = squeeze(mean(freq_avg.powspctrm, 1));
end

tim_interp = linspace(cfg_main.toi(1), cfg_main.toi(2), 512);
freq_interp = linspace(1, cfg_main.fmax, 512);

% We need to make a full time/frequency grid of both the original and
% interpolated coordinates. Matlab's meshgrid() does this for us:
[tim_grid_orig, freq_grid_orig] = meshgrid(tfr.time, tfr.freq);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

% And interpolate:
pow_interp = interp2(tim_grid_orig, freq_grid_orig, meanpow, tim_grid_interp, freq_grid_interp, 'spline');

%%
% while n==1
pow_interp1  = pow_interp(50:end,50:end);
tim_interp1  = tim_interp(50:end);
freq_interp1 = freq_interp(50:end);

%%
if cfg_main.effect == 2
    [~,idx] = min(pow_interp1(:));
else
    [~,idx] = max(pow_interp1(:));
end
[row,col] = ind2sub(size(pow_interp1),idx);

time_of_interest = tim_interp1(col);
freq_of_interest = freq_interp1(row);

timind = nearest(tim_interp, time_of_interest);
freqind = nearest(freq_interp, freq_of_interest);
pow_at_toi = pow_interp(:,timind);
pow_at_foi = pow_interp(freqind,:);

%%
figure();
ax_main  = axes('Position', [0.1 0.2 0.55 0.55]);
ax_right = axes('Position', [0.7 0.2 0.1 0.55]);
ax_top   = axes('Position', [0.1 0.8 0.55 0.1]);

%%
axes(ax_main);
im_main = imagesc(tim_interp, freq_interp, pow_interp);
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

%%
axes(ax_top);
area(tim_interp, pow_at_foi,...
    'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
xlim([cfg_main.toi(1), cfg_main.toi(2)]);
ylim([-clim clim]);
box off;
ax_top.XTickLabel = [];
ylabel('Power (dB)');
hold on;
plot([0 0], [-clim clim], 'k:');

%%
axes(ax_right);
area(freq_interp, pow_at_toi,...
    'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
view([270 90]); % this rotates the plot
ax_right.YDir = 'reverse';
ylim([-clim clim]);
box off;
ax_right.XTickLabel = [];
ylabel('Power (dB)');

%%
h = colorbar(ax_main, 'manual', 'Position', [0.85 0.2 0.05 0.55]);
ylabel(h, 'Power vs baseline (dB)');

%%
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

%%
cfg = [];
cfg.baseline = [-0.3 0];
cfg.baselinetype = 'absolute';
freq_bsl = ft_freqbaseline(cfg, tfr);

cfg = [];
cfg.variance = 'yes';
freq_sem = ft_freqdescriptives(cfg, freq_bsl);

tscore = freq_sem.powspctrm./ freq_sem.powspctrmsem;

% Average the t-score over our channels:
tscore = squeeze(mean(tscore, 1));

tscore(isnan(tscore))=0;

tscore_interp = interp2(tim_grid_orig, freq_grid_orig, tscore,...
    tim_grid_interp, freq_grid_interp, 'spline');

alpha = 0.01;
tcrit = tinv(1-alpha/2, size(tfr.powspctrm, 1)-1);

opacity = abs(tscore_interp) / tcrit;
opacity(opacity > 1) = 1;

if ~isempty(cfg_main.savefile)
    % hcp_write_figure([cfg_main.savefile,'.png'], gcf, 'resolution', 300);
    saveas(gcf,[cfg_main.savefile,'.png'])
end
end

function ep_data = vy_epoch(r_data,toi)

ep_data.all = r_data;

cfg = [];
cfg.toilim = toi(1,:);
ep_data.bsl = ft_redefinetrial(cfg, r_data);

cfg.toilim = toi(2,:);
ep_data.pst = ft_redefinetrial(cfg, r_data);

end

function [freq,ff, psd,tapsmofrq] = vy_fft(cfg_mian, data)

%%
cfg              = [];
cfg.method       = 'mtmfft';
cfg.output       = 'fourier';
cfg.keeptrials   = 'yes';
cfg.foilim       = cfg_mian.foilim;
cfg.tapsmofrq    = cfg_mian.tapsmofrq;
cfg.taper        = cfg_mian.taper; %'hanning';
% cfg.taper        = 'dpss';
cfg.pad          = 4;
freq             = ft_freqanalysis(cfg, data);
psd = squeeze(mean(mean(abs(freq.fourierspctrm),2),1));
ff = linspace(1, cfg.foilim(2), length(psd));

if cfg_mian.plotflag ==1
    figure,plot(ff,psd)
    xlabel('Hz'); ylabel('psd')
end

tapsmofrq = cfg.tapsmofrq;

if cfg_mian.saveflag ==1
    hcp_write_figure([cfg_mian.savefile,'.png'], gcf, 'resolution', 300);
end
end

function stat = vy_source_stat_montcarlo(s_data)

cfg = [];
cfg.parameter        = 'pow';
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'fdr';
cfg.clusteralpha     = 0.001;
% cfg.clusterstatistic = 'maxsum';
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