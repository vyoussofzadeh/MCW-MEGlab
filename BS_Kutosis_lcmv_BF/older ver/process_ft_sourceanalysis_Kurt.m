function varargout = process_ft_sourceanalysis_Kurt(varargin )
% PROCESS_FT_SOURCEANALYSIS Call FieldTrip function ft_sourceanalysis

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2019 University of Southern California & McGill University
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
% Authors: Francois Tadel, 2016-2017
% Vahab Youssof Zadeh, 2020-2021 % -- Adding DICS beamformer.

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% ===== PROCESS =====
% Description the process
sProcess.Comment     = 'FieldTrip: ft_sourceanalysis kurtosis v063021';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 356;
sProcess.Description = 'https://github.com/vyoussofzadeh/DICS-beamformer-for-Brainstorm';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data'};
sProcess.OutputTypes = {'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Label: Warning
sProcess.options.label1.Comment = '<B>Warning</B>: this is test, process under development.<BR><BR>';
sProcess.options.label1.Type    = 'label';
% Option: Inverse method
sProcess.options.method.Comment = 'Inverse method:';
sProcess.options.method.Type    = 'combobox_label';
sProcess.options.method.Value   = {'mne', {'Kurtosis LCMV beamformer', 'SAM beamformer', 'DICS beamformer', 'MNE', 'sLORETA', 'eLORETA', 'MUSIC', 'PCC', 'Residual variance'; ...
    'lcmv',            'sam',            'dics',            'mne', 'sloreta', 'eloreta', 'music', 'pcc', 'rv'}};
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
% Initialize fieldtrip
bst_ft_init();

h = get(0,'Children');
close(h)

%%
% ===== GET OPTIONS =====
% Inverse options
Method   = sProcess.options.method.Value{1};
Modality = sProcess.options.sensortype.Value{1};
% Get unique channel files
AllChannelFiles = unique({sInputs.ChannelFile});
% Progress bar
bst_progress('start', 'ft_sourceanalysis', 'Loading input files...', 0, 2*length(sInputs));

% ===== LOOP ON FOLDERS =====
for iChanFile = 1:1%length(AllChannelFiles)
    bst_progress('text', 'Loading input files...');
    % Get the study
    [sStudyChan, ~] = bst_get('ChannelFile', AllChannelFiles{iChanFile});
    % Error if there is no head model available
    if isempty(sStudyChan.iHeadModel)
        bst_report('Error', sProcess, [], ['No head model available in folder: ' bst_fileparts(sStudyChan.FileName)]);
        continue;
    end
    % Load channel file
    ChannelMat = in_bst_channel(AllChannelFiles{iChanFile});
    % Get selected sensors
    iChannels = channel_find(ChannelMat.Channel, Modality);
    if isempty(iChannels)
        bst_report('Error', sProcess, sInput, ['Channels "' Modality '" not found in channel file.']);
        return;
    end
    % Load head model
    HeadModelFile = sStudyChan.HeadModel(sStudyChan.iHeadModel).FileName;
    HeadModelMat = in_bst_headmodel(HeadModelFile);
    
    % ===== LOOP ON DATA FILES =====
    % Get data files for this channel file
    iChanInputs = find(ismember({sInputs.ChannelFile}, AllChannelFiles{iChanFile}));
    % Loop on data files
    for iInput = 1:1%length(iChanInputs)
        
        % === LOAD DATA ===
        % Load data
        DataFile = sInputs(iChanInputs(iInput)).FileName;
        DataMat = in_bst_data(DataFile);
        iStudyData = sInputs(iChanInputs(iInput)).iStudy;
        % Remove bad channels
        iBadChan = find(DataMat.ChannelFlag == -1);
        iChannelsData = setdiff(iChannels, iBadChan);
        % Error: All channels tagged as bad
        if isempty(iChannelsData)
            bst_report('Error', sProcess, sInput, 'All the selected channels are tagged as bad.');
            return;
        end
        % Convert data file to FieldTrip format
        ftData = out_fieldtrip_data(DataMat, ChannelMat, iChannelsData, 1);
        % Add data covariance
        %         ftData.cov = NoiseCovMat.NoiseCov(iChannelsData,iChannelsData);
        % Convert head model to FieldTrip format
        [ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);
        
        % === CALLING FIELDTRIP FUNCTION ===
        bst_progress('text', 'Calling FieldTrip function: ft_sourceanalysis...');
        % Prepare FieldTrip cfg structure
        cfg           = [];
        cfg.method    = Method;
        cfg.grid      = ftLeadfield;
        cfg.headmodel = ftHeadmodel;
        % Additional options for the method
        
        %% Step1, initial settings
        iChanFile = 1;
        ChannelMat = in_bst_channel(AllChannelFiles{iChanFile});
        iChannels = channel_find(ChannelMat.Channel, Modality);
        iChanInputs = find(ismember({sInputs.ChannelFile}, AllChannelFiles{iChanFile}));
        
        %% Step2, reading trials
        for iInput = 1:length(iChanInputs)
            DataFile = sInputs(iChanInputs(iInput)).FileName;
            DataMat = in_bst_data(DataFile);
            iStudyData = sInputs(iChanInputs(iInput)).iStudy;
            iBadChan = find(DataMat.ChannelFlag == -1);
            iChannelsData = setdiff(iChannels, iBadChan);
            
            if isempty(iChannelsData)
                bst_report('Error', sProcess, sInput, 'All the selected channels are tagged as bad.');
                return;
            end
            trl{iInput} = DataMat.F(iChannelsData,:);
            timee {iInput} = DataMat.Time;
        end
        
        ftData = out_fieldtrip_data(DataMat, ChannelMat, iChannelsData, 1);
        ftData1 = [];
        ftData1.label = ftData.label;
        ftData1.dimord = ftData.dimord;
        ftData1.trial = trl;
        ftData1.time = timee;
        
        %% step3, headmodel & leadfields ..
        [sStudyChan, ~] = bst_get('ChannelFile', AllChannelFiles{iChanFile});
        HeadModelFile = sStudyChan.HeadModel(sStudyChan.iHeadModel).FileName;
        HeadModelMat = in_bst_headmodel(HeadModelFile);
        Index = strfind(HeadModelMat.SurfaceFile, '/');
        subj = HeadModelMat.SurfaceFile(1:Index(1)-1);
        
        
        OutputDir = bst_fileparts(file_fullpath(DataFile));
        Index = strfind(OutputDir, 'data/');
        bsdir = OutputDir(1:Index(end)-1);
        bsdatadir = fullfile(bsdir,'data');
        bsanatdir = fullfile(bsdir,'anat');
        cd(bsdir)
        
        %% step4, saving path
        Index = strfind(DataFile, '/');
        saveid = DataFile(Index(end)+1:end-4);
        savepath = fullfile([subj, '_dics_',saveid]);
        if exist(savepath, 'file') == 0, mkdir(savepath), end
        
        %% Preprocessing
        datain = ftData1;
        
        disp('hpfreq to lpfreq eg, [18,25] in Hz:');
        foi = input('');
        cfg        = [];
        cfg.demean = 'yes';
        cfg.dftfilter = 'yes';
        cfg.hpfilter = 'yes';
        cfg.lpfilter = 'yes';
        cfg.hpfiltord = 3;
        cfg.lpfreq = foi(2);
        cfg.hpfreq = foi(1);
        fcln_data        = ft_preprocessing(cfg, datain);
        
        %%
        cfg = [];
        cfg.resamplefs = 500;
        data_resampled = ft_resampledata(cfg, fcln_data);
        
        %% ICA cleaning
        cfg            = [];
        cfg.method     = 'runica';
        cfg.numcomponent = 20;       % specify the component(s) that should be plotted
        comp           = ft_componentanalysis(cfg, data_resampled);
        
        %%
        cfg = [];
        cfg.layout = 'neuromag306mag.lay';
        lay = ft_prepare_layout(cfg);
        
        cfg = [];
        cfg.viewmode = 'component';
        cfg.layout = lay;
        ft_databrowser(cfg, comp);
        
        %%
        disp('Select bad ICs for:');
        bic = input('');
        close all;
        
        if ~ isempty(bic)
            
            cfg = [];
            cfg.component = comp.label(bic);
            cfg.updatesens = 'no';
            cln_data = ft_rejectcomponent(cfg, comp, data_resampled);
            
        else
            disp('no correction was done')
            cln_data = data_resampled;
        end
        
        %% Cutting 10 sec from the end of data, to reduce artifacts (eg., button press).
%         ct = 1;
%         if ct==1
%             cfg = [];
%             cfg.toilim = [cln_data.time{:}(10*cln_data.fsample),cln_data.time{:}(end-10*cln_data.fsample)];
%             cln_data = ft_redefinetrial(cfg,cln_data);
%         end
        
        %%
        cfg = [];
        cfg.channel = 'MEG';
        cfg.covariance = 'yes';
        cov_matrix = ft_timelockanalysis(cfg, cln_data);
        cov_matrix.grad = ftData.grad;
        
        %% Step8: head (forward) model
        sourcemodel = ft_read_headshape(fullfile(bsanatdir,HeadModelMat.SurfaceFile));
        [ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);
        
        %% step9: Surface-based source analysis
        cfg = [];
        cfg.method = 'lcmv';
        cfg.sourcemodel = ftLeadfield;
        cfg.headmodel = ftHeadmodel;
        cfg.lcmv.keepfilter = 'yes';
        cfg.lcmv.fixedori = 'yes'; % project on axis of most variance using SVD
        cfg.lcmv.lambda = '5%';
        %                  cfg.lcmv.lambda = '100%';
        cfg.lcmv.kappa = 69;
        cfg.lcmv.projectmom = 'yes'; % project dipole time series in direction of maximal power (see below)
        cfg.lcmv.kurtosis = 'yes';
        source = ft_sourceanalysis(cfg, cov_matrix);
        
        array = source.avg.kurtosis;
        array(isnan(array)) = 0;
        ispeak = imregionalmax(array); % findpeaksn is an alternative that does not require the image toolbox
        peakindex = find(ispeak(:));
        [~, i] = sort(source.avg.kurtosis(peakindex), 'descend'); % sort on the basis of kurtosis value
        peakindex = peakindex(i);
        
        npeaks = 3;
        disp(source.pos(peakindex(1:npeaks),:)); % output the positions of the top peaks
        
        
        %% Marking potential spikes in the source time series
        dat = ft_fetch_data(cln_data);
        hdr = ft_fetch_header(cln_data);
        
        for i=1:size(dat,1)
            %             hdr.label{i}= ['S' num2str(i)];
            hdr.chantype{i} = 'MEG';
            hdr.chanunit{i} = 'T' ; % see note below about scaling
        end
        
        % npeaks = 5;
        for i = 1:npeaks
            dat(end+1,:) = source.avg.mom{peakindex(i),:}; % see comment below about scaling
            hdr.label{end+1}= ['S' num2str(i)];
            hdr.chantype{end+1} = 'Source';
            hdr.chanunit{end+1} = 'T' ; % see note below about scaling
        end
        hdr.nChans = hdr.nChans+npeaks;
        % ft_write_data('Case3_timeseries', dat, 'header', hdr, 'dataformat', 'anywave_ades');
        dat1 = dat;
        
        k = 1;
        time_occur = [];
        dat = source.avg.mom{peakindex(1),:};
        sd = std(dat);
        PKS = [];
        tr = zeros(size(dat));
        kk = 20;
        while length(PKS) < 3
            tr(dat>kk*sd)=1;
            [PKS, peaksample] = findpeaks(tr, 'MinPeakDistance', 300); % peaks have to be separated by 300 sample points to be treated as separate
            kk = kk-1;
        end
        for j = 1:length(peaksample)
            time_occur(j) = source.time(peaksample(j));
            k = k + 1;
        end
        % disp(time_occur)
        
        toccur = [];
        for j=1:length(time_occur)
            toccur{j} = [num2str(j), '=', num2str(time_occur(j))];
        end
        disp(toccur')
        
        %%
        kk = 1;
        un_time_occur = unique(time_occur);
        L = length(un_time_occur);
        dat(peaksample)
        disp(toccur')
        keep_matrix = 1:length(toccur);
        un_time_occur_keep = un_time_occur(keep_matrix);
        
        %%
        mask = 'kurtosis';
        %- Spikes source time series plotting
        data_dummy = [];
        data_dummy.trial{1} = dat1;
        data_dummy.time = cln_data.time;
        data_dummy.hdr = hdr;
        data_dummy.label = hdr.label;
        
        clc
        tsel = [];
        for j=1:length(un_time_occur_keep)
            tsel{j} = [num2str(j), '=', num2str(un_time_occur_keep(j))];
        end
        disp(tsel')
        
        ETS_ask = 1;
        
        %     tselin = input('choose timing:');
        for j=1:length(un_time_occur_keep)
            disp([num2str(j), '/', num2str(length(un_time_occur_keep))])
            %
            cfg = [];
            cfg.toilim = [un_time_occur_keep(j) - kk,un_time_occur_keep(j) + kk];
            data_spk = ft_redefinetrial(cfg, cln_data);
            
            cfg = [];
            cfg.channel = 'MEG';
            cfg.covariance = 'yes';
            cov_matrix1 = ft_timelockanalysis(cfg, data_spk);
            cov_matrix1.grad = ftData.grad;
            
            cfg = [];
            cfg.method = 'lcmv';
            cfg.sourcemodel = ftLeadfield;
            cfg.headmodel = ftHeadmodel;
            cfg.lcmv.keepfilter = 'yes';
            cfg.lcmv.fixedori = 'yes'; % project on axis of most variance using SVD
            cfg.lcmv.lambda = '5%';
            %                  cfg.lcmv.lambda = '100%';
            cfg.lcmv.kappa = 69;
            cfg.lcmv.projectmom = 'yes'; % project dipole time series in direction of maximal power (see below)
            cfg.lcmv.kurtosis = 'yes';
            source_sel = ft_sourceanalysis(cfg, cov_matrix1);
            s_all(j,:) = source_sel.avg.(mask);
        end
        
        % === CREATE OUTPUT STRUCTURE ===
        bst_progress('text', 'Saving source file...');
        bst_progress('inc', 1);
        % Create structure
        ResultsMat = db_template('resultsmat');
        ResultsMat.ImagingKernel = [];
        switch Method
            case 'lcmv'
                ResultsMat.ImageGridAmp  = mean(s_all,1)';
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
        ResultsMat.Comment       = ['ft_sourceanalysis_kurtlcmv: ' Method, '_', datestr(now, 'dd/mm/yy-HH:MM')];
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
        ResultFile = bst_process('GetNewFilename', OutputDir, ['results_', Method, '_', Modality, ]);
        % Save new file structure
        bst_save(ResultFile, ResultsMat, 'v6');
        
        % ===== REGISTER NEW FILE =====
        bst_progress('inc', 1);
        % Create new results structure
        newResult = db_template('results');
        newResult.Comment       = ResultsMat.Comment;
        newResult.FileName      = file_short(ResultFile);
        newResult.DataFile      = DataFile;
        newResult.isLink        = 0;
        newResult.HeadModelType = ResultsMat.HeadModelType;
        % Get output study
        sStudyData = bst_get('Study', iStudyData);
        % Add new entry to the database
        iResult = length(sStudyData.Result) + 1;
        sStudyData.Result(iResult) = newResult;
        % Update Brainstorm database
        bst_set('Study', iStudyData, sStudyData);
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