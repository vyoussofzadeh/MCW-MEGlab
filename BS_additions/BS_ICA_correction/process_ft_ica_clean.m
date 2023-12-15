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

% Option: Number of ICA components
sProcess.options.icanum.Comment = 'Number of ICA components:';
sProcess.options.icanum.Type    = 'value';
sProcess.options.icanum.Value   = {20, 'components', 0}; % Default value is 20

% Option: Sensors selection
sProcess.options.lay.Comment = 'layout:';
sProcess.options.lay.Type    = 'combobox_label';
sProcess.options.lay.Value   = {'neuromag', {'neuromag', '4D', 'ctf', 'nolay'; ...
    'neuromag','4D','CTF', 'no layout'}};

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
% Using the user-defined ICA number or default if not provided
cfg.numcomponent = sProcess.options.icanum.Value{1};
cfg.method     = 'runica';
comp           = ft_componentanalysis(cfg, ftData1);

%% LAYOUT
layout = sProcess.options.lay.Value{1};

cfg = [];
switch layout
    case 'neuromag'
        cfg.layout = 'neuromag306mag.lay';
        lay = ft_prepare_layout(cfg);
    case '4D'
        cfg.layout = '4D248.lay';
        lay = ft_prepare_layout(cfg);
    case 'cft'
        disp('TO DO ..')
end

cfg = [];
switch layout
    case {'neuromag', '4D'}
        cfg.layout = lay;
end

cfg.viewmode = 'component';
% do_ft_databrowser(cfg, comp);
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
n = sProcess.options.icanum.Value{1};

% Calculate subplot dimensions
nby1 = ceil(sqrt(n));
nby2 = ceil(n / nby1);


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
close all;

if ~ isempty(bic)
    cfg = [];
    cfg.component = comp.label(bic);
    cfg.updatesens = 'no';
    ic_data = ft_rejectcomponent(cfg, comp, ftData1);
    
    %%
    ProtocolInfo = bst_get('ProtocolInfo');
    [datapath,~] = fileparts(DataFile);
    iStudyOut = sInputs(1).iStudy;
    
    OutputFiles = [];
    cd(fullfile(ProtocolInfo.STUDIES, datapath))
    for IInput = 1:length(sInputs)
        DataFile = sInputs(IInput).FileName;
        D = load(file_fullpath(DataFile));
        [A,FileName] = fileparts(DataFile);
        FileName_new = [FileName, '_ic'];
        tkz = tokenize(D.Comment, ' ');
        D.Comment = [tkz{1} '_ica ',tkz{2}];
        D.F(iChannelsData,:) = ic_data.trial{IInput};
        save(FileName_new, '-struct', 'D');
        
        newData = db_template('data');
        newData.Comment       = [tkz{1} '_ica ',tkz{2}];
        newData.FileName      = fullfile(A, [FileName_new, '.mat']);
        % Get output study
        sStudyOut = bst_get('Study', iStudyOut);
        % Add new entry to the database
        iResult = length(sStudyOut.Data) + 1;
        sStudyOut.Data(iResult) = newData;
        % Update Brainstorm database
        bst_set('Study', iStudyOut, sStudyOut);
        % Store output filename
        OutputFiles{end+1} = newData.FileName;
    end
    disp('done, reload datafile!')
    
else
    disp('no correction was applied')
    OutputFiles = [];
    for IInput = 1:length(sInputs)
        OutputFiles{end+1} = sInputs(IInput).FileName;
    end
end

db_save();
bst_progress('stop');

end
