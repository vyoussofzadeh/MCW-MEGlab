function varargout = process_ft_sourceanalysis_conn(varargin )
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
% Authors: Vahab YoussofZadeh, 2022

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'FieldTrip: ft_connanalysis (cross-spectral density)';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 357;
sProcess.Description = 'https://github.com/vyoussofzadeh/DICS-beamformer-for-Brainstorm';
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

    % Label: Frequency
    sProcess.options.label2.Comment = '<BR><B>Conn resolution:</B>';
    sProcess.options.label2.Type    = 'label';
    % Enter the FOI in the data in Hz, eg, 22:
    sProcess.options.conn.Comment = 'Conn res:';
    sProcess.options.conn.Type    = 'value';
    sProcess.options.conn.Value   = {1500, 'voxels (sruf points)', 0};

% Effects
sProcess.options.fconn.Comment = 'full resolution';
sProcess.options.fconn.Type    = 'checkbox';
sProcess.options.fconn.Value   = 1;

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
Connres = sProcess.options.conn.Value{1};
fconn = sProcess.options.fconn.Value;

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

%%
cfg = [];
cfg.toilim = PostStim;
ep_data = ft_redefinetrial(cfg, ftData);

%%
t_data = do_timelock(ep_data);

%%
cfg                  = [];
cfg.method           = 'lcmv';
cfg.sourcemodel  = ftLeadfield;
cfg.headmodel = ftHeadmodel;
cfg.keepfilter       = 'yes';
cfg.lcmv.keepfilter  = 'yes';
cfg.keeptrials       = 'yes';
cfg.lcmv.fixedori    = 'yes'; % project on axis of most variance using SVD
cfg.lcmv.lambda      = '0.1%';
source = ft_sourceanalysis(cfg, t_data);

if fconn == 1
    Connres  =length(source.avg.mom);
end

idx = round(linspace(1,length(source.avg.mom),Connres));
% idx = round(linspace(1,length(source.avg.mom),length(source.avg.mom)));
mom = []; k=1;
for i=idx
    clc
    disp([num2str(i),'/', num2str(length(source.avg.mom))])
    mom(k,:) = source.avg.mom{i};
    k = k+1;
end

outsum = do_conn(mom);
v = eigenvector_centrality_und(outsum);
% figure,bar(v), title('eigenvector');

ImageGridAmp = zeros(length(source.avg.mom),1);
for i=1:length(idx)-1
    ImageGridAmp(idx(i):idx(i+1)) = v(i);    
end

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


Method = 'conn';
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
% ResultsMat.Comment       = 'Conn_PLV';
ResultsMat.Comment       = ['ft_connanalysis: CrossSpectral_' Method, '_', datestr(now, 'dd/mm/yy-HH:MM')];
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

%
cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
cfg.preproc.demean   = 'yes';    % enable demean to remove mean value from each single trial
t_data            = ft_timelockanalysis(cfg, data);

end

function [outsum] = do_conn(mom)

[~, nrpt] = size(mom);
crsspctrm = (mom*mom')./nrpt;
tmp = crsspctrm; crsspctrm = []; crsspctrm(1,:,:) = tmp;

input = crsspctrm;
pownorm = 1;

siz = [size(input) 1];
% crossterms are described by chan_chan_therest
outsum = zeros(siz(2:end));
outssq = zeros(siz(2:end));
% outcnt = zeros(siz(2:end));
for j = 1:siz(1)
    if pownorm
        p1  = zeros([siz(2) 1 siz(4:end)]);
        p2  = zeros([1 siz(3) siz(4:end)]);
        for k = 1:siz(2)
            p1(k,1,:,:,:,:) = input(j,k,k,:,:,:,:);
            p2(1,k,:,:,:,:) = input(j,k,k,:,:,:,:);
        end
        p1    = p1(:,ones(1,siz(3)),:,:,:,:);
        p2    = p2(ones(1,siz(2)),:,:,:,:,:);
        denom = sqrt(p1.*p2); clear p1 p2;
    end
    tmp    = abs(reshape(input(j,:,:,:,:,:,:), siz(2:end))./denom); % added this for nan support marvin
    %tmp(isnan(tmp)) = 0; % added for nan support
    outsum = outsum + tmp;
    outssq = outssq + tmp.^2;
%     outcnt = outcnt + double(~isnan(tmp));
end

% size(outsum)
figure,imagesc(outsum), colorbar, title('conn (across voxels)');

end

function   v = eigenvector_centrality_und(CIJ)
%EIGENVECTOR_CENTRALITY_UND      Spectral measure of centrality
%
%   v = eigenvector_centrality_und(CIJ)
%
%   Eigenector centrality is a self-referential measure of centrality:
%   nodes have high eigenvector centrality if they connect to other nodes
%   that have high eigenvector centrality. The eigenvector centrality of
%   node i is equivalent to the ith element in the eigenvector 
%   corresponding to the largest eigenvalue of the adjacency matrix.
%
%   Inputs:     CIJ,        binary/weighted undirected adjacency matrix.
%
%   Outputs:      v,        eigenvector associated with the largest
%                           eigenvalue of the adjacency matrix CIJ.
%
%   Reference: Newman, MEJ (2002). The mathematics of networks.
%
%   Contributors:
%   Xi-Nian Zuo, Chinese Academy of Sciences, 2010
%   Rick Betzel, Indiana University, 2012
%   Mika Rubinov, University of Cambridge, 2015

%   MODIFICATION HISTORY
%   2010/2012: original (XNZ, RB)
%   2015: ensure the use of leading eigenvector (MR)


n = length(CIJ);
if n < 1000
    [V,D] = eig(CIJ);
else
    [V,D] = eigs(sparse(CIJ));
end
[~,idx] = max(diag(D));
ec = abs(V(:,idx));
v = reshape(ec, length(ec), 1);

end

function [ftHeadmodel, ftLeadfield, iChannels] = out_fieldtrip_headmodel_edt(HeadModelFile, ChannelFile, iChannels, isIncludeRef)
% OUT_FIELDTRIP_HEADMODEL: Converts a head model file into a FieldTrip structure (see ft_datatype_headmodel).
% 
% USAGE:  [ftHeadmodel, ftLeadfield, iChannels] = out_fieldtrip_headmodel(HeadModelFile, ChannelFile, isIncludeRef=1);
%         [ftHeadmodel, ftLeadfield, iChannels] = out_fieldtrip_headmodel(HeadModelMat,  ChannelMat,  isIncludeRef=1);
%
% INPUTS:
%    - HeadModelFile  : Relative path to a head model file available in the database
%    - HeadModelMat   : Brainstorm head model file structure
%    - ChannelFile    : Relative path to a channel file available in the database
%    - ChannelMat     : Brainstorm channel file structure
% OUTPUTS:
%    - ftHeadmodel    : Volume conductor model, typically returned by ft_prepare_headmodel
%    - ftLeadfield    : Leadfield matrix, typically returned by ft_prepare_leadfield
%    - iChannels      : Modified list of channels (after adding channels)

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
% Authors: Jeremy T. Moreau, Elizabeth Bock, Francois Tadel, 2015

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
    
    %% This section was diasbled, 05/11/22
%     iRef = channel_find(ChannelMat.Channel, 'MEG REF');
%     if ~isempty(iRef)
%         iChannels = [iRef, iChannels];
%     end
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
