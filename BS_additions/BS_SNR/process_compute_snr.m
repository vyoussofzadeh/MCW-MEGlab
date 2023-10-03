function varargout = process_compute_snr(varargin )
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
sProcess.Comment     = 'Computer SNR';
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
Modality = sProcess.options.sensortype.Value{1};

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
if isempty(iChannelsData)
    bst_report('Error', sProcess, sInputs, 'All the selected channels are tagged as bad.');
    return;
elseif any(~isChannelGood)
    bst_report('Info', sProcess, sInputs, ['Found ' num2str(length(find(~isChannelGood))) ' bad channels: ', sprintf('%s ', ChannelMat.Channel(find(~isChannelGood)).Name)]);
end

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

%% Calculate SNR

for j=1:size(ftData_A.trial,2)
    D1 = ftData_B.trial{j};
    D2 = ftData_A.trial{j};
    for i=1:size(D2,1)
        corrected_data = D2(i,:);
        original_data = D1(i,:);
        
        % Compute variances
        signal_power = var(corrected_data(:));
        noise_power = var(original_data(:) - corrected_data(:));
        
        % Calculate SNR
        SNR(j,i) = 10 * log10(signal_power / noise_power);
    end
end

mSNR = mean(mean(SNR,2));

figure,
subplot 221
imagesc(SNR), xlabel(['Sensor ', Modality]), ylabel({'Trial'}), title(['mean SNR (dB):', num2str(mSNR)])
colorbar
subplot 222
plot(mean(SNR,1),'.'), xlabel('Sensor'), ylabel({'mean SNR (dB)'}), set(gca,'color','none'), title(['min:', num2str(min(mean(SNR,1))), ', max:', num2str(max(mean(SNR,1)))])
subplot 223
imagesc(SNR'), xlabel('Trial'), ylabel({'Sensor'}), title(['mean SNR (dB):', num2str(mSNR)])
colorbar
subplot 224
plot(mean(SNR,2),'.'), xlabel('Trial'), ylabel({'mean SNR (dB)'}),set(gca,'color','none'), title(['min:', num2str(min(mean(SNR,2))), ', max:', num2str(max(mean(SNR,2)))])
set(gcf, 'Position', [600   300   600   600]);


% Create a new data structure for SNR
sDataOut = db_template('data');
sDataOut.F            = SNR';  % Transpose SNR so that each column corresponds to a trial
sDataOut.Comment      = 'SNR';
sDataOut.ChannelFlag  = ones(length(ChannelMat.Channel), 1);
sDataOut.Time         = 1:size(SNR, 1);  % Each entry corresponds to a trial
sDataOut.DataType     = 'recordings';    % This is continuous data
sDataOut.Device       = 'SNR';
sDataOut.nAvg         = 1;               % Number of averaged epochs
sDataOut.DisplayUnits = 'dB';

% Get output study
[sStudy, iStudy] = bst_get('AnyFile', sInputsA(1).FileName);

% Save the new file in Brainstorm database
OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_snr');
sDataOut.FileName = file_short(OutputFile);
bst_save(OutputFile, sDataOut, 'v7');

% Update database
db_add_data(iStudy, OutputFile, sDataOut);

% % Save database
db_save();
% % Hide progress bar
bst_progress('stop');

end


