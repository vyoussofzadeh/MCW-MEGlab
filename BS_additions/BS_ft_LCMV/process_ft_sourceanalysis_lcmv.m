function varargout = process_ft_sourceanalysis_lcmv(varargin )
% process_ft_sourceanalysis_conn Call FieldTrip function ft_sourceanalysis (LCMV)

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
% Author: Vahab YoussofZadeh, 2022

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
% Description the process
sProcess.Comment     = 'FieldTrip: ft_sourceanalysis (LCMV)';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 357;
sProcess.Description = 'https://github.com/vyoussofzadeh/DICS-beamformer-for-Brainstorm';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data'};
sProcess.OutputTypes = {'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

sProcess.options.label1.Comment = '<BR><B>Time of interest:</B>';
sProcess.options.label1.Type    = 'label';
% Active time window
sProcess.options.poststim.Comment = 'Active (post-stim):';
sProcess.options.poststim.Type    = 'poststim';
sProcess.options.poststim.Value   = [];
% Baseline time window
sProcess.options.baseline.Comment = 'Baseline (pre-stim):';
sProcess.options.baseline.Type    = 'baseline';
sProcess.options.baseline.Value   = [];

% Option: Sensors selection
sProcess.options.sensortype.Comment = 'Sensor type:';
sProcess.options.sensortype.Type    = 'combobox_label';
sProcess.options.sensortype.Value   = {'MEG', {'MEG', 'MEG GRAD', 'MEG MAG', 'EEG', 'SEEG', 'ECOG'; ...
    'MEG', 'MEG GRAD', 'MEG MAG', 'EEG', 'SEEG', 'ECOG'}};

sProcess.options.erds.Comment = {'ERD', 'ERS', 'both', 'ERD/S effects:'; ...
    'erd', 'ers', 'both', ''};
sProcess.options.erds.Type    = 'radio_linelabel';
sProcess.options.erds.Value   = 'erd';
% Effects
sProcess.options.effect.Comment = {'abs', 'raw', 'Absolute/raw value of the contrast:'; ...
    'abs', 'raw', ''};
sProcess.options.effect.Type    = 'radio_linelabel';
sProcess.options.effect.Value   = 'abs';

% % Label: Frequency
% sProcess.options.label3.Comment = '<BR><B>Conn resolution:</B>';
% sProcess.options.label3.Type    = 'label';
% % Enter the FOI in the data in Hz, eg, 22:
% sProcess.options.conn.Comment = 'Conn res:';
% sProcess.options.conn.Type    = 'value';
% sProcess.options.conn.Value   = {1500, 'voxels (sruf points)', 0};
%
% % Effects
% sProcess.options.fconn.Comment = 'full resolution';
% sProcess.options.fconn.Type    = 'checkbox';
% sProcess.options.fconn.Value   = 1;

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
Baseline = sProcess.options.baseline.Value{1};
Modality = sProcess.options.sensortype.Value{1};
% Connres = sProcess.options.conn.Value{1};
% fconn = sProcess.options.fconn.Value;

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
[ftHeadmodel, ftLeadfield, iChannelsData] = out_fieldtrip_headmodel_edt(HeadModelMat, ChannelMat, iChannelsData, 1);

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

%% A workaround to avoid using ft_chanunit
unit = []; for i=1:length(ftData.grad.chantype); unit{i} = 'meg'; end
ftData.grad.chanunit = unit';

%%
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
    ftData        = ft_preprocessing(cfg, ftData);
else
    ftData = ftData;
end

%%
% Baseline
cfg = []; cfg.toilim = Baseline;
ep_data.bsl = ft_redefinetrial(cfg, ftData);
% Post-stim
cfg.toilim = PostStim;
ep_data.pst = ft_redefinetrial(cfg, ftData);
% Baseline + Post-stim
cfg = [];
ep_data.app = ft_appenddata(cfg, ep_data.bsl, ep_data.pst);

%%
cov_matrix = do_timelock(ep_data);

%%
cfg = [];
cfg.method = 'lcmv';
cfg.sourcemodel  = ftLeadfield;
cfg.headmodel = ftHeadmodel;
cfg.lcmv.lambda = '100%';
% cfg.lcmv.lambda     = '0.1%';
cfg.lcmv.keepfilter = 'yes';
source_whole = ft_sourceanalysis(cfg, cov_matrix.app);

cfg =[];
cfg.method = 'lcmv';
cfg.sourcemodel  = ftLeadfield;
cfg.sourcemodel.filter = source_whole.avg.filter;
cfg.headmodel = ftHeadmodel;
% cfg.rawtrial = 'yes';
% cfg.keeptrials = 'yes';
% cfg.lcmv.projectmom = 'yes';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.lambda     = '10%';
s_data.bsl      = ft_sourceanalysis(cfg, cov_matrix.bsl);
s_data.pst      = ft_sourceanalysis(cfg, cov_matrix.pst);

% s_data.pst.pow = s_data.pst.avg.pow;
% s_data.bsl.pow = s_data.bsl.avg.pow;

% s_data.pst.avg.pow - s_data.bsl.avg.pow

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
%%
% ===== SAVE RESULTS =====
% === CREATE OUTPUT STRUCTURE ===
bst_progress('text', 'Saving source file...');
bst_progress('inc', 1);
% Output study
if (length(sInputs) == 1)
    iStudyOut = sInputs(1).iStudy;
    RefDataFile = sInputs(iChanInputs(iInput)).FileName;
else
    [~, iStudyOut] = bst_process('GetOutputStudy', sProcess, sInputs);
    RefDataFile = [];
end
% Create structure
ResultsMat = db_template('resultsmat');
ResultsMat.ImagingKernel = [];

ImageGridAmp = source_diff_dics.pow;

Method = 'lcmv';
ResultsMat.ImageGridAmp  = ImageGridAmp;
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
ResultsMat.Comment       = ['ft_sourceanalysis:' Method, '_', num2str(foi(1)),'_',num2str(foi(2)),'_Hz_', datestr(now, 'dd/mm/yy-HH:MM')];
switch lower(ResultsMat.HeadModelType)
    case 'volume'
        ResultsMat.GridLoc    = HeadModelMat.GridLoc;
    case 'surface'
        ResultsMat.GridLoc    = [];
    case 'mixed'
        ResultsMat.GridLoc    = HeadModelMat.GridLoc;
        ResultsMat.GridOrient = HeadModelMat.GridOrient;
end
ResultsMat = bst_history('add', ResultsMat, 'compute', ['ft_connanalysis: ' Method ' ' Modality ' ']);

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

function t_data = do_timelock(data)

cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
cfg.preproc.demean   = 'yes';    % enable demean to remove mean value from each single trial
cfg.keeptrials       = 'yes';
t_data.pst            = ft_timelockanalysis(cfg, data.pst);
t_data.bsl            = ft_timelockanalysis(cfg, data.bsl);
cfg.vartrllength     = 2;
t_data.app           = ft_timelockanalysis(cfg, data.app);

end


function [ftHeadmodel, ftLeadfield, iChannels] = out_fieldtrip_headmodel_edt(HeadModelFile, ChannelFile, iChannels, isIncludeRef)
% ===== PARSE INPUTS =====
if (nargin < 4) || isempty(isIncludeRef)
    isIncludeRef = 1;
end
ftHeadmodel = [];
ftLeadfield = [];

% ===== LOAD INPUTS =====
% Load head model file
if ischar(HeadModelFile)
    HeadModelMat = in_bst_headmodel(HeadModelFile);
elseif isstruct(HeadModelFile)
    HeadModelMat = HeadModelFile;
else
    error('Failed to load head model.');
end
% If this file was computed with FieldTrip, it should include the original FieldTrip headmodel
if isfield(HeadModelMat, 'ftHeadmodelMeg') && ~isempty(HeadModelMat.ftHeadmodelMeg)
    ftHeadmodel = HeadModelMat.ftHeadmodelMeg;
elseif isfield(HeadModelMat, 'ftHeadmodelEeg') && ~isempty(HeadModelMat.ftHeadmodelEeg)
    ftHeadmodel = HeadModelMat.ftHeadmodelEeg;
end
% Load channel file
if ischar(ChannelFile)
    ChannelMat = in_bst_channel(ChannelFile);
elseif isstruct(ChannelFile)
    ChannelMat = ChannelFile;
else
    error('Failed to load channel file.')
end
% Get sensor type
Modality = ChannelMat.Channel(iChannels(1)).Type;
% Get headmodel type
switch (Modality)
    case 'EEG',   HeadModelMethod = HeadModelMat.EEGMethod;
    case 'ECOG',  HeadModelMethod = HeadModelMat.ECOGMethod;
    case 'SEEG',  HeadModelMethod = HeadModelMat.SEEGMethod;
    case {'MEG','MEG MAG','MEG GRAD','MEG REF'}, HeadModelMethod = HeadModelMat.MEGMethod;
end

% ===== ADD MEG REF =====
if isIncludeRef && ismember(Modality, {'MEG','MEG MAG','MEG GRAD'})
end


% ===== CREATE FIELDTRIP HEADMODEL =====
if isempty(ftHeadmodel)
    % Get subject
    sSubject = bst_get('SurfaceFile', HeadModelMat.SurfaceFile);
    % Headmodel type
    switch (HeadModelMethod)
        case {'meg_sphere', 'singlesphere'}
            ftHeadmodel.type = 'singlesphere';
            ftHeadmodel.r = HeadModelMat.Param(iChannels(1)).Radii(1);
            ftHeadmodel.o = HeadModelMat.Param(iChannels(1)).Center(:)';
            
        case {'eeg_3sphereberg', 'concentricspheres'}
            ftHeadmodel.type = 'concentricspheres';
            ftHeadmodel.r = HeadModelMat.Param(iChannels(1)).Radii(:)';
            ftHeadmodel.o = HeadModelMat.Param(iChannels(1)).Center(:)';
            % Get default conductivities
            BFSProperties = bst_get('BFSProperties');
            ftHeadmodel.c = BFSProperties(1:3);
            
        case {'os_meg', 'localspheres'}
            ftHeadmodel.type = 'localspheres';
            ftHeadmodel.r = [HeadModelMat.Param(iChannels).Radii]';
            ftHeadmodel.o = [HeadModelMat.Param(iChannels).Center]';
            
        case {'singleshell'}
            ftHeadmodel.type = HeadModelMethod;
            % Check if the surfaces are available
            if isempty(sSubject.iInnerSkull)
                error('No inner skull surface available for this subject.');
            else
                disp(['BST> ' HeadModelMethod ': Using the default inner skull surface available in the database.']);
            end
            % Load surfaces
            SurfaceFiles = {sSubject.Surface(sSubject.iInnerSkull).FileName};
            ftHeadmodel.bnd = out_fieldtrip_tess(SurfaceFiles);
            
        case {'openmeeg', 'dipoli', 'bemcp'}
            ftHeadmodel.type = HeadModelMethod;
            % Check if the surfaces are available
            if isempty(sSubject.iInnerSkull) || isempty(sSubject.iOuterSkull) || isempty(sSubject.iScalp)
                error('No BEM surfaces available for this subject.');
            else
                disp(['BST> ' HeadModelMethod ': Using the default surfaces available in the database (inner skull, outer skull, scalp).']);
            end
            % Load surfaces
            SurfaceFiles = {sSubject.Surface(sSubject.iScalp).FileName, ...
                sSubject.Surface(sSubject.iOuterSkull).FileName, ...
                sSubject.Surface(sSubject.iInnerSkull).FileName};
            ftHeadmodel.bnd = out_fieldtrip_tess(SurfaceFiles);
            % Default OpenMEEG options
            ftHeadmodel.cond         = [0.33, 0.004125, 0.33];
            ftHeadmodel.skin_surface = 1;
            ftHeadmodel.source       = 3;
            % ERROR: MISSING .mat field??
            
        otherwise
            error('out_fieldtrip_headmodel does not support converting this type of head model.');
    end
    % Unit and labels
    ftHeadmodel.unit  = 'm';
    ftHeadmodel.label = {ChannelMat.Channel(iChannels).Name}';
end

% ===== CREATE FIELDTRIP LEADFIELD =====
if (nargout >= 2)
    % Create FieldTrip structure
    nSources = length(HeadModelMat.GridLoc);
    ftLeadfield.pos             = HeadModelMat.GridLoc;
    ftLeadfield.unit            = 'm';
    ftLeadfield.inside          = true(nSources, 1);
    ftLeadfield.leadfielddimord = '{pos}_chan_ori';
    ftLeadfield.label           = {ChannelMat.Channel(iChannels).Name};
    ftLeadfield.leadfield       = cell(1, nSources);
    for i = 1:nSources
        ftLeadfield.leadfield{i} = HeadModelMat.Gain(iChannels, 3*(i-1)+[1 2 3]);
    end
end
end
