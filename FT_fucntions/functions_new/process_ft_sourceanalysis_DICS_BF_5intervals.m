
function varargout = process_ft_sourceanalysis_DICS_BF_5intervals(varargin )
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
sProcess.Comment     = 'FieldTrip: ft_sourceanalysis tv-DICS-BF EML (earlylate) v031121';
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
sProcess.options.method.Value   = {'mne', {'LCMV beamformer', 'SAM beamformer', 'DICS beamformer', 'MNE', 'sLORETA', 'eLORETA', 'MUSIC', 'PCC', 'Residual variance'; ...
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
bst_progress('stop');

% ===== GET OPTIONS =====
% Inverse options
if iscell(sProcess.options.method.Value)
    Method   = sProcess.options.method.Value{1};
    Modality = sProcess.options.sensortype.Value{1};
else
    Method   = sProcess.options.method.Value;
    Modality = sProcess.options.sensortype.Value;
end
% Get unique channel files
AllChannelFiles = unique({sInputs.ChannelFile});
% Progress bar
% bst_progress('start', 'ft_sourceanalysis', 'Loading input files...', 0, 2*length(sInputs));

% ===== LOOP ON FOLDERS =====
for iChanFile = 1:1%length(AllChannelFiles)
    bst_progress('text', 'Loading input files...');
    % Get the study
    [sStudyChan, ~] = bst_get('ChannelFile', AllChannelFiles{iChanFile});
    % Error if there is no head model available
    if isempty(sStudyChan.iHeadModel)
        bst_report('Error', sProcess, [], ['No head model available in folder: ' bst_fileparts(sStudyChan.FileName)]);
        continue;
    elseif isempty(sStudyChan.NoiseCov) || isempty(sStudyChan.NoiseCov(1).FileName)
        bst_report('Error', sProcess, [], ['No noise covariance matrix available in folder: ' bst_fileparts(sStudyChan.FileName)]);
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
    % Load data covariance matrix
    NoiseCovFile = sStudyChan.NoiseCov(1).FileName;
    NoiseCovMat = load(file_fullpath(NoiseCovFile));
    %%% DATA OR NOISE COVARIANCE ????
    
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
        ftData.cov = NoiseCovMat.NoiseCov(iChannelsData,iChannelsData);
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
        % Loop on data files
        for iInput = 1:length(iChanInputs)
            % === LOAD DATA ===
            % Load data
            DataFile = sInputs(iChanInputs(iInput)).FileName;
            DataMat = in_bst_data(DataFile);
            iStudyData = sInputs(iChanInputs(iInput)).iStudy;
            % Remove bad channels
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
        ftData1.grad = ftData.grad;
        ftData1.dimord = ftData.dimord;
        ftData1.trial = trl;
        ftData1.time = timee;
        
        %% step3, headmodel & leadfields ..
        [sStudyChan, ~] = bst_get('ChannelFile', AllChannelFiles{iChanFile});
        % Load head model
        HeadModelFile = sStudyChan.HeadModel(sStudyChan.iHeadModel).FileName;
        HeadModelMat = in_bst_headmodel(HeadModelFile);
        % Load data covariance matrix
        NoiseCovFile = sStudyChan.NoiseCov(1).FileName;
        NoiseCovMat = load(file_fullpath(NoiseCovFile));
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
        savepath = fullfile(['dics_',saveid]);
%         if exist(savepath, 'file') == 0, mkdir(savepath), end
        
        %% step 5: Selecting freq of interest
        if exist(fullfile(savepath,'tfr.mat'),'file')~= 2
            
            %             disp('Enter the highest freq in the data in Hz, eg, 40:')
            %             fmax = input(':');
            fmax = 40;
            % do tfr-decomposition
            cfg = [];
            cfg.output     = 'pow';
            cfg.channel    = 'all';
            cfg.method     = 'mtmconvol';
            cfg.taper      = 'hanning';
            % cfg.taper      = 'dpss';
            cfg.foi        = 1:2:fmax;
            cfg.keeptrials = 'yes';
            cfg.t_ftimwin  = 3./cfg.foi;
            % cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
            cfg.tapsmofrq  = 0.8 *cfg.foi;
            cfg.toi        = ftData1.time{1}(1):0.05:ftData1.time{1}(end);
            tfr            = ft_freqanalysis(cfg, ftData1);
%             save(fullfile(savepath,'tfr.mat'),'tfr','fmax');
            
%         else
%             load(fullfile(savepath,'tfr.mat'));
%             if ~exist('fmax','var'), fmax = 30; end
        end
        
        %         cfg = [];
        %         cfg.savepath = 1;
        %         cfg.savefile = fullfile(savepath,'tfr');
        %         cfg.fmax = fmax;
        %         cfg.toi = [ftData1.time{1}(1), ftData1.time{1}(end)];
        %         [time_of_interest,freq_of_interest] = vy_tfr_plot(cfg, tfr);
        %         disp(['time_of_interest:',num2str(time_of_interest),'sec']);
        %         disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
%         L = 0.3;
        
        %% Time varying window
        w1 = 0; l = 0.8; ov = l.*0.2; j=1; wi=[];
%         w1 = 0.4; l = 0.8; ov = l.*0.2; j=1; wi=[];
        while w1+l < 2
            wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
        end
        
        %% Step 6: selecting time of interest
        datain = ftData1;

        %% Step8: head (forward) model
        sourcemodel = ft_read_headshape(fullfile(bsanatdir,HeadModelMat.SurfaceFile));
        [ftHeadmodel, ftLeadfield] = out_fieldtrip_headmodel(HeadModelMat, ChannelMat, iChannelsData, 1);
        %%
        dics = [];
%         f = 22; tapsmofrq = 4;
        f = 20; tapsmofrq = 5;
        
        %% Selecting trials from 3 time intervals
        tkz = tokenize(OutputDir,'/'); tkz = tokenize(tkz{end},'_');
        outdir = '/data/MEG/Research/ECP/';
        outd.sub = fullfile(outdir,'/Story_Math/MEG_work_ft',tkz{2}, tkz{1}, tkz{3});
        
        load(fullfile(outd.sub,['ic_',tkz{2}]));
        
        idx_str_sent1 = cln_data.trialinfo(:,21)==10;
        tmp  = cln_data.trialinfo(idx_str_sent1,8);
        idx_str_sent1_narration = unique(tmp);
        
%         tmp1 = tmp;
%         for i=1:length(idx_str_sent1_narration)
%             idx = tmp1 == idx_str_sent1_narration(i);
%             tmp1(idx) = i;
%         end
        
        tmp1 = tmp;
        for i=1:length(idx_str_sent1_narration)
            idx = find(tmp == idx_str_sent1_narration(i));
            tmp1(idx) = i;
        end
        
        clear l
        for j=1:length(idx_str_sent1_narration)
            l(j) = floor(length(find(tmp1 == j))/4);
        end
        
        sel = [];
        sel.first = []; sel.seond = []; sel.third = []; sel.fourth = []; sel.fifth = [];
        for k=1:length(l)
            if isempty(find(l(k) == 0, 1))
%                 disp(l(k))
%                 pause,
                idx = find(tmp1 == k);
                sel.first = [sel.first, idx(1):idx(l(k))];
                sel.seond = [sel.seond, idx(l(k))+1:idx(l(k))+l(k)];
                sel.third = [sel.third, idx(l(k))+2:idx(l(k))+l(k)+1];
                sel.fourth = [sel.fourth, idx(l(k))+3:idx(l(k))+l(k)+2];
                sel.fifth = [sel.fifth, idx(l(k))+4:idx(end)]; %idx(l(k))+l(k)+3 
            end
        end
        
        %%  Trial selection
        t_sel = sProcess.options.t_sel.Value;
        switch t_sel
            case 1
                trl_sel = sel.first; slabel = '1st';
            case 2
                trl_sel = sel.seond; slabel = '2nd';
            case 3
                trl_sel = sel.third; slabel = '3rd';
            case 4
                trl_sel = sel.fourth; slabel = '4th';
            case 5
                trl_sel = sel.fifth; slabel = '5th';
        end
        
        for i=1:size(wi,1)
            
            toi(1,:) = [-0.3,0];
            toi(2,:) = wi(i,:);
            ep_data = vy_epoch(datain, toi);
            ep_data_all  = ep_data.pst;
            
            ep_data1 = [];
            cfg = [];
            cfg.trials = trl_sel;
            ep_data1.pst = ft_redefinetrial(cfg, ep_data.pst);
            ep_data1.bsl = ft_redefinetrial(cfg, ep_data.bsl);
            
            cfg = [];
            %             cfg.trials = trl_sel(1,:);
            ep_data1.bsl = ep_data1.bsl;
            ep_data = ep_data1;
            ep_data.all = ep_data_all;
            
            %%
            cfg = [];
            ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.all);
            
            %% step7: spectral analysis
            cfg_main = [];
            %             cfg_main.fmax = fmax;
            cfg_main.sens = datain.grad;
            cfg_main.outputdir = savepath;
            %             cfg_main.freq_of_interest  = freq_of_interest; % Hz
            
            cfg = [];
            cfg.savefile = [];
            cfg.saveflag = 2;
            cfg.foilim = [f f];
            cfg.plotflag  = 2;
            cfg.taper    = 'dpss'; cfg.tapsmofrq  = tapsmofrq;
            
            if f < 4, cfg.tapsmofrq  = 1; cfg.taper    = 'hanning'; end
            
            [f_data.app,~,~,~] = vy_fft(cfg, ep_data.app); f_data.app.elec = cfg_main.sens;
            f_data.bsl = vy_fft(cfg, ep_data.bsl); f_data.bsl.elec = cfg_main.sens;
            f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = cfg_main.sens;
            
            %%
            %             disp('Simple post-pre contrasting (fast)         : 1');
            %             disp('Cluster-based permutation statistics (slow): 2');
            %             disp('?');
            %             st = input('');
            st = 1;
            if st == 2, Method = 'dics_stat'; end
            
            %% Step9: Source analysis
            switch Method
                case 'mne'
                    cfg.mne.prewhiten = 'yes';
                    cfg.mne.lambda    = 3;
                    cfg.mne.scalesourcecov = 'yes';
                    Time = DataMat.Time;
                    
                case 'lcmv'
                    Time = [DataMat.Time(1), DataMat.Time(2)];
                    
                case 'dics'
                    %- FFT_based
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
                    %                 disp('desynchronisation effects: 1, synchronisation effects: 2:');
                    %                 ask_sd = input(':');
                    ask_sd = 1;
                    
                    switch ask_sd
                        case 1
                            %
                            cfg = [];
                            cfg.parameter = 'pow';
                            cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
                            source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
                            source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
                            %                             source_diff_dics.pow(source_diff_dics.pow>0)=0;
                            %                             source_diff_dics.pow = abs(source_diff_dics.pow);
                            
                        case 2
                            cfg = [];
                            cfg.parameter = 'pow';
                            cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
                            source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
                            source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
                            %                             source_diff_dics.pow(source_diff_dics.pow<0)=0;
                            %                             source_diff_dics.pow = abs(source_diff_dics.pow);
                            
                            %
                            %                         figure
                            %                         m = source_diff_dics.pow;
                            %                         bnd.pnt = sourcemodel.pos;
                            %                         bnd.tri = sourcemodel.tri;
                            %                         ft_plot_mesh(bnd, 'vertexcolor', abs(m));
                            %                         colorbar
                    end
                case 'dics_stat'
                    
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
                    
                    stat = vy_source_stat_montcarlo(s_data);
                    
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
            switch Method
                case 'dics'
                    dics{i} = source_diff_dics;
                case 'dics_stat'
                    stats_all{i} = stats2;
            end
        end
        
        %% Plotting mean
        switch Method
            case 'dics'
                for i=1:length(dics)
                    dicspow(i,:) =  dics{i}.pow;
                end
                dicspow(isnan(dicspow))=0;
                
                if size(dicspow,1) > 1
                    dicsmean = mean(dicspow);
                else
                    dicsmean = dicspow;
                end
                
            case 'dics_stat'
                for i=1:length(stats_all)
                    stat_value(i,:) =  stats_all{i}.stat;
                end
                stat_value(isnan(stat_value))=0;
                
                if size(stat_value,1) > 1
                    stat_value = mean(stat_value);
                else
                    stat_value = stat_value;
                end
                
        end
        
        %% step10: saving ourput, surface/volume projection
        % Call FieldTrip function
        %             ftSource = ft_sourceanalysis(cfg, ftData);
        
        % === CREATE OUTPUT STRUCTURE ===
        bst_progress('text', 'Saving source file...');
        bst_progress('inc', 1);
        % Create structure
        ResultsMat = db_template('resultsmat');
        ResultsMat.ImagingKernel = [];
        switch Method
            case 'dics'
                ResultsMat.ImageGridAmp  = dicsmean';
                ResultsMat.cfg           = source_diff_dics.cfg;
            case 'dics_stat'
                ResultsMat.ImageGridAmp  = stat_value';
                ResultsMat.cfg           = stat.cfg;
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
        ResultsMat.Comment       = ['tv_EML_2sided_ft_sourceanalysis: ',slabel,'_', Method, '_',num2str(f),'Hz_', datestr(now, 'dd/mm/yy-HH:MM')];
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

% And baseline-correct the average:
cfg = [];
cfg.baseline = [-0.3 0];
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

%%
% while n==1
pow_interp1  = pow_interp(50:end,50:end);
tim_interp1  = tim_interp(50:end);
freq_interp1 = freq_interp(50:end);

%%
[~,idx] = min(pow_interp1(:));
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