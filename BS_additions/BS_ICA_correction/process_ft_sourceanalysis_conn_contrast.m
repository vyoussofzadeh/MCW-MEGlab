function varargout = process_ft_sourceanalysis_conn_contrast(varargin )
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
sProcess.Comment     = 'FieldTrip: ft_connanalysis contrast (cross-spectral density)';
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
% Baseline time window
sProcess.options.baseline.Comment = 'Baseline (pre-stim):';
sProcess.options.baseline.Type    = 'baseline';
sProcess.options.baseline.Value   = [];

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
Baseline = sProcess.options.baseline.Value{1};
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
ep = []; cfg = [];
cfg.toilim = PostStim;
ep.pst = ft_redefinetrial(cfg, ftData);
cfg.toilim = Baseline;
ep.bsl = ft_redefinetrial(cfg, ftData);

cfg = [];
ep.app = ft_appenddata(cfg,ep.bsl,ep.pst);

%%
% tlk = [];
% tlk.pst = do_timelock(ep.pst);
% tlk.bsl = do_timelock(ep.bsl);
% tlk.app = do_timelock(ep.app);

%%
cfg                  = [];
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
cfg.preproc.demean   = 'yes';    % enable demean to remove mean value from each single trial
cfg.keeptrials       = 'yes';
% cfg.vartrllength = 0;

tlk.pst            = ft_timelockanalysis(cfg, ep.pst);
tlk.bsl            = ft_timelockanalysis(cfg, ep.bsl);

cfg.vartrllength     = 2;
tlk.app           = ft_timelockanalysis(cfg, ep.app);

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
source = ft_sourceanalysis(cfg, tlk.app);

cfg = [];
cfg.method           = 'lcmv';
cfg.sourcemodel        = ftLeadfield;
cfg.sourcemodel.filter = source.avg.filter;
cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
cfg.headmodel = ftHeadmodel;
src.bsl      = ft_sourceanalysis(cfg, tlk.bsl);
src.pst      = ft_sourceanalysis(cfg, tlk.pst);

if fconn == 1
    Connres  = length(src.pst.avg.mom);
end

idx = round(linspace(1,length(src.pst.avg.mom),Connres));
% idx = round(linspace(1,length(source.avg.mom),length(source.avg.mom)));

%%
mom = []; k=1;
for i=idx, clc
    disp([num2str(i),'/', num2str(length(src.pst.avg.mom))]); mom(k,:) = src.pst.avg.mom{i}; k = k+1;
end
conn.pst = do_conn(mom);

mom = []; k=1;
for i=idx, clc, disp([num2str(i),'/', num2str(length(src.bsl.avg.mom))])
    mom(k,:) = src.bsl.avg.mom{i}; k = k+1;
end
conn.bsl = do_conn(mom);

%%
net.pst = eigenvector_centrality_und(conn.pst);
net.bsl = eigenvector_centrality_und(conn.bsl);
net_diff = net.pst - net.bsl;

% figure,bar(net_diff), title('eigenvector');

%%
ImageGridAmp = zeros(length(source.avg.mom),1);
for i=1:length(idx)-1
    ImageGridAmp(idx(i):idx(i+1)) = net_diff(i);    
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
    [tmp, iStudyOut] = bst_process('GetOutputStudy', sProcess, sInputs);
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
% figure,imagesc(outsum), colorbar, title('conn (across voxels)');

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
