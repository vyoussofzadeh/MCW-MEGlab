function varargout = process_computeLI(varargin )
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
% Authors: Vahab YoussofZadeh, 2023

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>


% Description the process
sProcess.Comment     = 'Surface-based (Desikan-Killiany) laterality';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 337;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/CoregisterSubjects';
sProcess.InputTypes  = {'results'};
sProcess.OutputTypes = {'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

% Add an option to select the LI computation method
sProcess.options.LImethod.Comment = 'Select LI computation method:';
sProcess.options.LImethod.Type    = 'combobox';
sProcess.options.LImethod.Value   = {1, {'Counting', 'Bootstrapping'}};

% Add a time interval input
sProcess.options.time_interval.Comment = 'Choose a time interval:';
sProcess.options.time_interval.Type    = 'combobox';
sProcess.options.time_interval.Value   = {1, {'Time interval', 'averaged sources'}};

% Add an effect input
sProcess.options.effect.Comment = 'Choose an effect:';
sProcess.options.effect.Type    = 'combobox';
sProcess.options.effect.Value   = {1, {'Positive values', 'Negative', 'Absolute'}};

% Add a threshold type input
sProcess.options.threshtype.Comment = 'Threshold type:';
sProcess.options.threshtype.Type    = 'combobox';
sProcess.options.threshtype.Value   = {1, {'Global-max: all time and regions combined', ...
    'Time-max: time of interest (toi) and all regions', ...
    'Region-max: toi and regions of interests (rois)'}};
% Add a threshold ratio input
sProcess.options.ratio4threshold.Comment = 'Threshold ratio (%):';
sProcess.options.ratio4threshold.Type    = 'value';
sProcess.options.ratio4threshold.Value   = {50, '%', 0, 100, 1, 1}; % {default value, '', min value, max value, step size, decimal digits}
% sProcess.options.ratio4threshold.Value   = {0.5};

% Add a folder directory input
sProcess.options.savedir.Comment = 'Saving Dir.:';
sProcess.options.savedir.Type    = 'text';
sProcess.options.savedir.Value   = ''; % Default value can be empty or a specific path

% Add a saving file name input
sProcess.options.sname.Comment = 'Saving filename:';
sProcess.options.sname.Type    = 'text';
sProcess.options.sname.Value   = ''; % Default value can be empty or a specific name

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput)

OutputFiles = {};
sResultP = in_bst_results(sInput.FileName, 1);

% Obtain saving directory
savedir = sProcess.options.savedir.Value;

% Prompt and select time interval
time_interval = selectTimeInterval(sProcess.options.time_interval.Value{1});

% Select data type and time range
timerange = selectDataTypeAndTimeRange(time_interval, sResultP);

% Select effect type (Positive, Negative, Absolute values)
effect = selectEffectType(sProcess.options.effect.Value{1});

% Determine sample rate
samplerate = determineSampleRate(time_interval, sResultP);

% Define ROI-related parameters
% TotROI = 8;  % Total number of ROIs for LI-calculation
% Ratio4Threshold = getThresholdRatio(sProcess.options.ratio4threshold.Value{1});
Ratio4Threshold = sProcess.options.ratio4threshold.Value{1}/100;

% Define threshold type
Threshtype = selectThresholdType(sProcess.options.threshtype.Value{1});

% Process ImageGridAmp based on selected effect
ImageGridAmp = processImageGridAmp(sResultP.ImageGridAmp, effect);

% Determine max values for various time windows
[AllMax, GlobalMax, t1, t2] = determineMaxValues(time_interval, ImageGridAmp, sResultP, timerange);

% Convert Desikan-Killiany scout to selected scouts
[sScout, ~] = convertDesikanKillianyScout(sResultP);

% Define ROIs
[RoiLabels, RoiIndices] = defineROIs();

% Load scout files and compute LI

cfg_LI = [];
cfg_LI.time_interval = time_interval;
cfg_LI.ImageGridAmp = ImageGridAmp;
cfg_LI.timerange = timerange;
cfg_LI.RoiLabels = RoiLabels;
cfg_LI.RoiIndices = RoiIndices;
cfg_LI.sScout = sScout;
cfg_LI.AllMax = AllMax;
cfg_LI.GlobalMax = GlobalMax;
cfg_LI.Threshtype = Threshtype;
cfg_LI.Ratio4Threshold = Ratio4Threshold;
cfg_LI.t1 = t1;
cfg_LI.t2 = t2;
cfg_LI.savedir = savedir;
cfg_LI.sname =  sProcess.options.sname.Value;

switch sProcess.options.LImethod.Value{1}
    case 1
        computeLI_counting(cfg_LI);
    case 2
        computeLI_bootstrap(cfg_LI);
end

disp('To edit the LI script, first ensure Brainstorm is running. Then, open process_computeLI.m in Matlab.');
disp('Pipeline update: 09/25/23');

end

% === HELPER FUNCTIONS ===
function time_interval = selectTimeInterval(time_interval)
% Prompt user to select time interval

% Ensure valid selection
if isempty(time_interval) || ~any(time_interval == [1, 2, 3])
    error('Invalid time interval selection. Choose 1, 2, or 3.');
end
end

function timerange = selectDataTypeAndTimeRange(time_interval, ~)
% Prompt user to select data type and determine the corresponding time range

switch time_interval
    case 1
        disp('Enter Timerange: e.g. [.100 .900]: ')
        timerange = input('');
        if isempty(timerange) || length(timerange) ~= 2
            error('Invalid Timerange input. Provide a valid range.');
        end
    case 2
        timerange = 1;  % Make sure these are less than 1 for mS
    otherwise
        error('Invalid time interval. Choose 1, 2, or 3.');
end
end

function effect = selectEffectType(effect)
% Prompt user to select effect type

% Ensure valid selection
if isempty(effect) || ~any(effect == [1, 2, 3])
    error('Invalid effect type selection. Choose 1, 2, or 3.');
end
end

function samplerate = determineSampleRate(time_interval, sResultP)
% Determine the sample rate based on the time interval and sResultP

switch time_interval
    case 2
        samplerate = 1;
    otherwise
        samplerate = round(inv((sResultP.Time(end) - sResultP.Time(1)) / length(sResultP.Time))) - 1;
end
end

function Threshtype = selectThresholdType(Threshtype)
% Prompt user to select threshold type

% Ensure valid selection
if isempty(Threshtype) || ~any(Threshtype == [1, 2, 3])
    error('Invalid threshold type selection. Choose 1, 2, or 3.');
end
end

function ImageGridAmp = processImageGridAmp(ImageGridAmp, effect)
% Apply the desired effect on the ImageGridAmp

switch effect
    case 1
        % Positive values: No change needed, as ImageGridAmp remains the same
    case 2
        ImageGridAmp = -ImageGridAmp;
    case 3
        ImageGridAmp = abs(ImageGridAmp);
    otherwise
        error('Invalid effect type. Choose 1 (Positive), 2 (Negative), or 3 (Absolute).');
end
end

function [AllMax, GlobalMax, t1, t2] = determineMaxValues(time_interval, ImageGridAmp, sResultP, timerange)
% Compute the maximum values AllMax and GlobalMax

switch time_interval
    case 2
        GlobalMax = max(ImageGridAmp(:));  % Max value over all time points
        AllMax = max(ImageGridAmp(:));     % Max value over the time window of interest
        t1 = []; t2 = [];
    otherwise
        t1 = find(sResultP.Time >= timerange(1), 1);
        t2 = find(sResultP.Time >= timerange(2), 1);
        AllMax = max(max(ImageGridAmp(:, t1:t2)));   % Max value over the time window of interest
        GlobalMax = max(max(ImageGridAmp));           % Max value over all time points
end
end

function [sScout, ProtocolInfo] = convertDesikanKillianyScout(sResultP)
% Convert the Desikan-Killiany scout to select scouts

ProtocolInfo = bst_get('ProtocolInfo');
SurfaceFile = load(fullfile(ProtocolInfo.SUBJECTS, sResultP.SurfaceFile));

Scouts = [];
sScout = [];
for i = 1:length(SurfaceFile.Atlas)
    if contains(SurfaceFile.Atlas(i).Name, {'Desikan-Killiany'})
        Scouts = SurfaceFile.Atlas(i).Scouts;
    end
end
sScout.Scouts = Scouts;

% Handle case when number of anatomical regions are not identical to atlas regions
l = length(sScout.Scouts);
if l ~= 68
    rois = {sScout.Scouts.Label};
    rois_temp = {sScout.Scouts(1:68).Label};
    [C, IA] = setdiff(rois_temp, rois);
    
    warning('The number of anatomical regions are not identical to atlas regions');
    disp('Replacing with zero ...');
    new_eScout = sScout;
    k = 1;
    for i = 1:68
        if i ~= IA
            new_eScout.Scouts(i) = sScout.Scouts(k);
            k = k + 1;
        else
            new_eScout.Scouts(i).Vertices = [];
            new_eScout.Scouts(i).Seed = [];
        end
    end
    sScout = new_eScout;
end
end

function [RoiLabels, RoiIndices] = defineROIs()
% Define regions of interest (ROIs)

AngSmg   = [15,16,63,64];
Front    = [3,4,5,6,11,12,25,26,29,30,33,34,37,38,39,40,41,42,49,50,53,54,55,56,57,58];
LatFront = [5,6,11,12,37,38,39,40,41,42,55,56,57,58];
LatTemp  = [1,2,17,18,31,32,61,62,65,66,67,68];
PeriSyl  = [15,16,37,38,41,42,61,62,63,64];
Tanaka   = [37,38,41,42,61,62,63,64];
Temp     = [1,2,9,10,13,14,17,18,19,20,27,28,31,32,35,36,61,62,65,66,67,68];
Whole    = 1:68;

RoiLabels = {'AngSmg', 'Front','LatFront','LatTemp', 'PeriSyl', 'Tanaka','Temp','Whole'};
RoiIndices = {AngSmg, Front, LatFront, LatTemp, PeriSyl, Tanaka, Temp, Whole};
end

function computeLI_counting(cfg_LI)
% Compute the Laterality Index (LI) and associated tasks

time_interval = cfg_LI.time_interval;
ImageGridAmp = cfg_LI.ImageGridAmp;
timerange = cfg_LI.timerange;
RoiLabels = cfg_LI.RoiLabels;
RoiIndices  = cfg_LI.RoiIndices;
sScout  = cfg_LI.sScout;
AllMax = cfg_LI.AllMax;
GlobalMax = cfg_LI.GlobalMax;
Threshtype = cfg_LI.Threshtype;
Ratio4Threshold = cfg_LI.Ratio4Threshold;
t1 = cfg_LI.t1;
t2 = cfg_LI.t2;
savedir = cfg_LI.savedir;

TotROI = 8;

s1='LI_';
Summ_LI=zeros(1,TotROI); % initialize the vector that summarizes the final LIs  % added JL 11212014
Summ_LI_Label='ROI Labels: '; % initialize the string that summarizes the ROI labels  % added JL 11212014
switch time_interval
    case 1
        figure
end
plot_ind=1;
LI_label_out={};
%
for ii = 1:8
    
    s2 = RoiLabels{ii};
    %Odd indices are left Rois
    Ltemp_region = [];
    Ltemp_label  = [];
    hemi_roi_num=length(RoiIndices{ii});
    curr_subregion=sScout.Scouts(RoiIndices{ii});
    
    k = 1;
    for i=1:2:hemi_roi_num
        Ltemp_region=[Ltemp_region,curr_subregion(i).Vertices];
        Ltemp_label{k} = curr_subregion(i).Label; k = k+1;
    end
    
    %Even indices are right Rois
    Rtemp_region = [];
    Rtemp_label  = [];
    k = 1;
    for i=2:2:hemi_roi_num
        Rtemp_region=[Rtemp_region,curr_subregion(i).Vertices];
        Rtemp_label{k} = curr_subregion(i).Label; k = k+1;
    end
    LHscout = Ltemp_region;
    RHscout = Rtemp_region;
    
    switch time_interval % modified by VY
        case 2
            %First parse the maps into separate space-times maps for each side
            LHvals = ImageGridAmp(LHscout);
            LH_max = max(max(LHvals));
            RHvals = ImageGridAmp(RHscout);
            RH_max = max(max(RHvals));
            ROIMax = max(LH_max,RH_max);
        otherwise
            %First parse the maps into separate space-times maps for each side
            LHvals = ImageGridAmp(LHscout,t1:t2);
            LH_max = max(max(LHvals));
            RHvals = ImageGridAmp(RHscout,t1:t2);
            RH_max = max(max(RHvals));
            ROIMax = max(LH_max,RH_max);
    end
    
    switch Threshtype %modified by vy@09/08/22
        case 1
            threshold = Ratio4Threshold*GlobalMax;  % dSPM threshold to get rid of non-significant voxels. Added JL@10/30/14
        case 2
            threshold = Ratio4Threshold*AllMax;  % dSPM threshold to get rid of non-significant voxels. Added JL@10/30/14
        case 3
            threshold = Ratio4Threshold*ROIMax; % Added by VY@09/08/22
    end
    
    ind_L = find(LHvals > threshold);
    ind_R = find(RHvals > threshold);
    
    % ROI count --- above threshold voxels only
    %JL@10/30/14    indHalfROImaxL = find(LHvals > ROIMax/2);%indHalfROImaxL = find(LHvals > ROIMax/2); % Not being used right now
    L_ROIcount = length(ind_L); L_count(ii) = L_ROIcount;
    %JL@10/30/14    indHalfROImaxR = find(RHvals > ROIMax/2); % Not being used right now
    R_ROIcount = length(ind_R); R_count(ii) = R_ROIcount;
    ROIcount=sum(L_ROIcount+R_ROIcount); % to report total significant voxels over space-time
    LI_ROIcount = 100*((L_ROIcount-R_ROIcount)/(L_ROIcount+R_ROIcount));
    Summ_LI(ii)=LI_ROIcount;  % added JL 11212014
    Summ_LI_Label=[Summ_LI_Label  sprintf('\t') s2];   % added JL 11212014
    LI_label_out=[LI_label_out, s2];
    
    % ROI average --- above threshold voxels only
    LHvals_aboveThreshold = LHvals(ind_L); % a 1-D matrix, no need for mean(mean()) later
    RHvals_aboveThreshold = RHvals(ind_R); % a 1-D matrix, no need for mean(mean()) later
    L_ROIavg=mean(LHvals_aboveThreshold);
    R_ROIavg=mean(RHvals_aboveThreshold);
    ROIavg = mean([L_ROIavg,R_ROIavg]);
    LI_ROIavg = 100*((L_ROIavg-R_ROIavg)/(L_ROIavg+R_ROIavg));
    
    % Run a loop to plot LIs based on space-time voxel count as a function of threshold
    k=0;
    Rng= threshold:0.2:AllMax; % Rng= threshold:1:AllMax; % modified JL@10/30/14
    for thrTmp = Rng
        ind = find(LHvals > thrTmp);
        L_ROIcount  = length(ind);
        ind = find(RHvals > thrTmp);
        R_ROIcount  = length(ind);
        k=k+1;
        Thrshd_LI_ROIcount(k,1)=thrTmp;
        if L_ROIcount+R_ROIcount~=0
            LI_ROIcount = 100*((L_ROIcount-R_ROIcount)/(L_ROIcount+R_ROIcount));
        else
            LI_ROIcount = inf;
        end
        Thrshd_LI_ROIcount(k,2)=LI_ROIcount;
    end
    switch time_interval
        case 1
            subplot(4,4,plot_ind)
            plot_ind=plot_ind+1;
            plot(Thrshd_LI_ROIcount(:,1),Thrshd_LI_ROIcount(:,2));
            title([RoiLabels{ii} ' Count-based']);
    end
    
    % Run a loop to plot LIs based on space-time voxel-magnitude average as a function of threshold. But this metrix does not seem
    % to be very useful because the LIs tend to be around 0.
    k=0;
    %JL@10/30/14   Rng= threshold:1:AllMax;
    for thrTmp = Rng
        ind = LHvals > thrTmp;
        LHvals_aboveThreshold = LHvals(ind); % a 1-D matrix, no need for mean(mean()) later
        ind = RHvals > thrTmp;
        RHvals_aboveThreshold = RHvals(ind); % a 1-D matrix, no need for mean(mean()) later
        if isempty(LHvals_aboveThreshold)
            L_ROIavg=0;  %to prevent error of "Warning: Divide by zero" when LHvals_aboveThreshold is an empty matrix and you are doing operation of mean(LHvals_aboveThreshold)
        else
            L_ROIavg=mean(LHvals_aboveThreshold);
        end
        if isempty(RHvals_aboveThreshold)
            R_ROIavg=0;  %to prevent error of "Warning: Divide by zero" when RHvals_aboveThreshold is an empty matrix and you are doing operation of mean(RHvals_aboveThreshold)
        else
            R_ROIavg=mean(RHvals_aboveThreshold);
        end
        
        k=k+1;
        Thrshd_LI_ROIavg(k,1)=thrTmp;
        if L_ROIavg+R_ROIavg ~= 0
            LI_ROIavg = 100*((L_ROIavg-R_ROIavg)/(L_ROIavg+R_ROIavg));
        else
            LI_ROIavg = inf;
        end
        Thrshd_LI_ROIavg(k,2)=LI_ROIavg;
    end
    
    switch time_interval
        case 1
            subplot(4,4,plot_ind)
            plot_ind=plot_ind+1;
            plot(Thrshd_LI_ROIavg(:,1),Thrshd_LI_ROIavg(:,2)); %JS 092815 changed plot to subplot
            title([RoiLabels{ii} ' Average-based']);
            set(gcf, 'Position', [500   500   1000   800]);
    end
end

% Save results to disk
% Create folder path if it doesn't exist
folderPath = fullfile(savedir);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
    disp('Folder created successfully.');
else
    disp('Folder already exists.');
end

sname = cfg_LI.sname;

% Define the filename based on the time_interval
if time_interval == 2
    filename = fullfile(folderPath, ['/LI_ROItable_', sname, '_thresh', num2str(Ratio4Threshold), '.xls']);
else
    filename = fullfile(folderPath, ['/LI_ROItable_', sname, '_', 'time', num2str(timerange(1)), '-', num2str(timerange(2)), '_thresh', num2str(Ratio4Threshold), '.xls']);
end

% Write the data to the Excel file
tempfile = fopen(filename, 'w');
fprintf(tempfile, '%s\t', LI_label_out{:});
fprintf(tempfile, '\n');
fprintf(tempfile, '%f\t', Summ_LI);
fprintf(tempfile, '\n\nThreshold');
fprintf(tempfile, '%f\t', threshold);
fclose(tempfile);

% Display the path to the saved file
disp(['Results saved to: ' filename]);

% Optional: Display results
disp('=================');
disp('                 ');
a = table(RoiLabels'); a.Properties.VariableNames{'Var1'} = 'ROI';
b = table(Summ_LI'); b.Properties.VariableNames{'Var1'} = 'LI';
d = [a, b];
disp(d);
disp('=================');

end


function computeLI_bootstrap(cfg_LI)
% Compute the Laterality Index (LI) and associated tasks

time_interval = cfg_LI.time_interval;
ImageGridAmp = cfg_LI.ImageGridAmp;
timerange = cfg_LI.timerange;
RoiLabels = cfg_LI.RoiLabels;
RoiIndices  = cfg_LI.RoiIndices;
sScout  = cfg_LI.sScout;
AllMax = cfg_LI.AllMax;
GlobalMax = cfg_LI.GlobalMax;
Threshtype = cfg_LI.Threshtype;
Ratio4Threshold = cfg_LI.Ratio4Threshold;
t1 = cfg_LI.t1;
t2 = cfg_LI.t2;
savedir = cfg_LI.savedir;

TotROI = 8;

s1='LI_';
Summ_LI=zeros(1,TotROI); % initialize the vector that summarizes the final LIs  % added JL 11212014
Summ_LI_Label='ROI Labels: '; % initialize the string that summarizes the ROI labels  % added JL 11212014
switch time_interval
    case 1
        figure
end
plot_ind=1;
LI_label_out={};
%
for ii = 1:8
    
    s2 = RoiLabels{ii};
    %Odd indices are left Rois
    Ltemp_region = [];
    Ltemp_label  = [];
    hemi_roi_num=length(RoiIndices{ii});
    curr_subregion=sScout.Scouts(RoiIndices{ii});
    
    k = 1;
    for i=1:2:hemi_roi_num
        Ltemp_region=[Ltemp_region,curr_subregion(i).Vertices];
        Ltemp_label{k} = curr_subregion(i).Label; k = k+1;
    end
    
    %Even indices are right Rois
    Rtemp_region = [];
    Rtemp_label  = [];
    k = 1;
    for i=2:2:hemi_roi_num
        Rtemp_region=[Rtemp_region,curr_subregion(i).Vertices];
        Rtemp_label{k} = curr_subregion(i).Label; k = k+1;
    end
    LHscout = Ltemp_region;
    RHscout = Rtemp_region;
    
    switch time_interval % modified by VY
        case 2
            %First parse the maps into separate space-times maps for each side
            LHvals = ImageGridAmp(LHscout);
            LH_max = max(max(LHvals));
            RHvals = ImageGridAmp(RHscout);
            RH_max = max(max(RHvals));
            ROIMax = max(LH_max,RH_max);
        otherwise
            %First parse the maps into separate space-times maps for each side
            LHvals = ImageGridAmp(LHscout,t1:t2);
            LH_max = max(max(LHvals));
            RHvals = ImageGridAmp(RHscout,t1:t2);
            RH_max = max(max(RHvals));
            ROIMax = max(LH_max,RH_max);
    end
    
%     %- Li Calc.
%     cfg = [];
%     cfg.sinput = S_data_sel.sFiles_in; %sFiles_anim_hc;
%     cfg.BS_data_dir = BS_data_dir;
%     cfg.atlas = Data_hcp_atlas.atlas; cfg.thre = thre; cfg.fplot = 0;
%     cfg.index_L = idx_L{j};
%     cfg.index_R = idx_R{j};
%     cfg.Threshtype = Threshtype;
%     cfg.divs = 100;
%     cfg.n_resampling = 200;
%     cfg.RESAMPLE_RATIO = 0.75;
%     cfg.dvd = 5; % 5
%     [weighted_li, LI_max, num_threshvals] = LI_bootstrap(cfg)

% end
divs = 100;

switch Threshtype %modified by vy@09/08/22
    case 1
        %             threshold = Ratio4Threshold*GlobalMax;  % dSPM threshold to get rid of non-significant voxels. Added JL@10/30/14
        threshold = (0:(divs-1)) * (GlobalMax / (divs - 1));
        
    case 2
        %             threshold = Ratio4Threshold*AllMax;  % dSPM threshold to get rid of non-significant voxels. Added JL@10/30/14
        threshold = (0:(divs-1)) * (AllMax / (divs - 1));
    case 3
        %             threshold = Ratio4Threshold*ROIMax; % Added by VY@09/08/22
        threshold = (0:(divs-1)) * (ROIMax / (divs - 1));
end


lr_threshvals = cell(1, numel(threshold));
num_threshvals = 0;

for i = 1:numel(threshold)
    threshval = threshvals(i);
    l_threshvals = lvals_nonnegative(lvals_nonnegative >= threshval);
    r_threshvals = rvals_nonnegative(rvals_nonnegative >= threshval);
    
    %- Check if both hemispheres have enough voxels above the threshold
    % (MIN_NUM_THRESH_VOXELS)
    if (numel(l_threshvals) >= MIN_NUM_THRESH_VOXELS && numel(r_threshvals) >= MIN_NUM_THRESH_VOXELS)
        num_threshvals = num_threshvals + 1;
        lr_threshvals{num_threshvals} = {l_threshvals, r_threshvals};
    else
        break;
    end
end
interval_length = length(lr_threshvals);

TB_LIs = zeros(num_threshvals, n_resampling);

for i = 1:num_threshvals
    l_set = lr_threshvals{i}{1};
    r_set = lr_threshvals{i}{2};
    
    l_n = max(round(numel(l_set) * RESAMPLE_RATIO), MIN_NUM_THRESH_VOXELS);
    r_n = max(round(numel(r_set) * RESAMPLE_RATIO), MIN_NUM_THRESH_VOXELS);
    
    l_indices = randi(numel(l_set), n_resampling, l_n);
    lactivity = sum(l_set(l_indices), 2);
    
    r_indices = randi(numel(r_set), n_resampling, r_n);
    ractivity = sum(r_set(r_indices), 2);
    
    TB_LIs(i, :) = (lactivity - ractivity) ./ (lactivity + ractivity);
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ind_L = find(LHvals > threshold);
    ind_R = find(RHvals > threshold);
    
    % ROI count --- above threshold voxels only
    %JL@10/30/14    indHalfROImaxL = find(LHvals > ROIMax/2);%indHalfROImaxL = find(LHvals > ROIMax/2); % Not being used right now
    L_ROIcount = length(ind_L); L_count(ii) = L_ROIcount;
    %JL@10/30/14    indHalfROImaxR = find(RHvals > ROIMax/2); % Not being used right now
    R_ROIcount = length(ind_R); R_count(ii) = R_ROIcount;
    ROIcount=sum(L_ROIcount+R_ROIcount); % to report total significant voxels over space-time
    LI_ROIcount = 100*((L_ROIcount-R_ROIcount)/(L_ROIcount+R_ROIcount));
    Summ_LI(ii)=LI_ROIcount;  % added JL 11212014
    Summ_LI_Label=[Summ_LI_Label  sprintf('\t') s2];   % added JL 11212014
    LI_label_out=[LI_label_out, s2];
    
    % ROI average --- above threshold voxels only
    LHvals_aboveThreshold = LHvals(ind_L); % a 1-D matrix, no need for mean(mean()) later
    RHvals_aboveThreshold = RHvals(ind_R); % a 1-D matrix, no need for mean(mean()) later
    L_ROIavg=mean(LHvals_aboveThreshold);
    R_ROIavg=mean(RHvals_aboveThreshold);
    ROIavg = mean([L_ROIavg,R_ROIavg]);
    LI_ROIavg = 100*((L_ROIavg-R_ROIavg)/(L_ROIavg+R_ROIavg));
    
    % Run a loop to plot LIs based on space-time voxel count as a function of threshold
    k=0;
    Rng= threshold:0.2:AllMax; % Rng= threshold:1:AllMax; % modified JL@10/30/14
    for thrTmp = Rng
        ind = find(LHvals > thrTmp);
        L_ROIcount  = length(ind);
        ind = find(RHvals > thrTmp);
        R_ROIcount  = length(ind);
        k=k+1;
        Thrshd_LI_ROIcount(k,1)=thrTmp;
        if L_ROIcount+R_ROIcount~=0
            LI_ROIcount = 100*((L_ROIcount-R_ROIcount)/(L_ROIcount+R_ROIcount));
        else
            LI_ROIcount = inf;
        end
        Thrshd_LI_ROIcount(k,2)=LI_ROIcount;
    end
    switch time_interval
        case 1
            subplot(4,4,plot_ind)
            plot_ind=plot_ind+1;
            plot(Thrshd_LI_ROIcount(:,1),Thrshd_LI_ROIcount(:,2));
            title([RoiLabels{ii} ' Count-based']);
    end
    
    % Run a loop to plot LIs based on space-time voxel-magnitude average as a function of threshold. But this metrix does not seem
    % to be very useful because the LIs tend to be around 0.
    k=0;
    %JL@10/30/14   Rng= threshold:1:AllMax;
    for thrTmp = Rng
        ind = LHvals > thrTmp;
        LHvals_aboveThreshold = LHvals(ind); % a 1-D matrix, no need for mean(mean()) later
        ind = RHvals > thrTmp;
        RHvals_aboveThreshold = RHvals(ind); % a 1-D matrix, no need for mean(mean()) later
        if isempty(LHvals_aboveThreshold)
            L_ROIavg=0;  %to prevent error of "Warning: Divide by zero" when LHvals_aboveThreshold is an empty matrix and you are doing operation of mean(LHvals_aboveThreshold)
        else
            L_ROIavg=mean(LHvals_aboveThreshold);
        end
        if isempty(RHvals_aboveThreshold)
            R_ROIavg=0;  %to prevent error of "Warning: Divide by zero" when RHvals_aboveThreshold is an empty matrix and you are doing operation of mean(RHvals_aboveThreshold)
        else
            R_ROIavg=mean(RHvals_aboveThreshold);
        end
        
        k=k+1;
        Thrshd_LI_ROIavg(k,1)=thrTmp;
        if L_ROIavg+R_ROIavg ~= 0
            LI_ROIavg = 100*((L_ROIavg-R_ROIavg)/(L_ROIavg+R_ROIavg));
        else
            LI_ROIavg = inf;
        end
        Thrshd_LI_ROIavg(k,2)=LI_ROIavg;
    end
    
    switch time_interval
        case 1
            subplot(4,4,plot_ind)
            plot_ind=plot_ind+1;
            plot(Thrshd_LI_ROIavg(:,1),Thrshd_LI_ROIavg(:,2)); %JS 092815 changed plot to subplot
            title([RoiLabels{ii} ' Average-based']);
            set(gcf, 'Position', [500   500   1000   800]);
    end
end

% Save results to disk
% Create folder path if it doesn't exist
folderPath = fullfile(savedir);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
    disp('Folder created successfully.');
else
    disp('Folder already exists.');
end

sname = cfg_LI.sname;

% Define the filename based on the time_interval
if time_interval == 2
    filename = fullfile(folderPath, ['/LI_ROItable_', sname, '_thresh', num2str(Ratio4Threshold), '.xls']);
else
    filename = fullfile(folderPath, ['/LI_ROItable_', sname, '_', 'time', num2str(timerange(1)), '-', num2str(timerange(2)), '_thresh', num2str(Ratio4Threshold), '.xls']);
end

% Write the data to the Excel file
tempfile = fopen(filename, 'w');
fprintf(tempfile, '%s\t', LI_label_out{:});
fprintf(tempfile, '\n');
fprintf(tempfile, '%f\t', Summ_LI);
fprintf(tempfile, '\n\nThreshold');
fprintf(tempfile, '%f\t', threshold);
fclose(tempfile);

% Display the path to the saved file
disp(['Results saved to: ' filename]);

% Optional: Display results
disp('=================');
disp('                 ');
a = table(RoiLabels'); a.Properties.VariableNames{'Var1'} = 'ROI';
b = table(Summ_LI'); b.Properties.VariableNames{'Var1'} = 'LI';
d = [a, b];
disp(d);
disp('=================');

end

function [weighted_li, LI_max, num_threshvals] = LI_bootstrap(cfg_main)

% MIN_NUM_THRESH_VOXELS is a predefined minimum number ...
% of voxels that are required for the computation to be considered valid.

wi = cfg_main.wi;
atlas = cfg_main.atlas;
idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;
divs = cfg_main.divs; % divs: number of divisions between 0 and max value in img
n_resampling = cfg_main.n_resampling; % 200

sinput = cfg_main.sinput;

RESAMPLE_RATIO = cfg_main.RESAMPLE_RATIO; % 0.75
MIN_NUM_THRESH_VOXELS = round(5 / RESAMPLE_RATIO);

tmp = load(fullfile(cfg_main.BS_data_dir, sinput));
sScout = atlas;

weighted_li = [];
% ft_progress('init', 'text',     'please wait ...');

for j = 1:size(wi, 1)
    
    %     ft_progress(j/size(wi, 1), 'Processing interval %d from %d', j, size(wi, 1));
    timind1 = nearest(tmp.Time, wi(j,1));
    timind2 = nearest(tmp.Time, wi(j,2));
    
    d_in = mean(tmp.ImageGridAmp(:, timind1:timind2), 2);
    ImageGridAmp = abs(d_in);
    
    % Get left and right subregions from scout data
    LHscout = [];
    for i = 1:length(idx_L)
        LHscout = [LHscout, sScout.Scouts(idx_L(i)).Vertices];
    end
    
    RHscout = [];
    for i = 1:length(idx_R)
        RHscout = [RHscout, sScout.Scouts(idx_R(i)).Vertices];
    end
    
    % Extract amplitude values for left and right subregions
    LHvals = ImageGridAmp(LHscout);
    RHvals = ImageGridAmp(RHscout);
    
    % Calculate maximum values for left and right subregions
    LH_max = max(LHvals(:));
    RH_max = max(RHvals(:));
    ROIMax = max(LH_max, RH_max);
    
    %     divs = 5; % divs: number of divisions between 0 and max value in img
    lvals_nonnegative = LHvals(LHvals >= 0);
    rvals_nonnegative = RHvals(RHvals >= 0);
    
    % creating an array of values between 0 and the maximum amplitude
    % (ROIMax) observed in the brain regions of interest.
    % This array is linearly spaced and divided into divs number of divisions.
    % This creates a series of thresholds that will be used to progressively ...
    % include only those voxels with amplitude values above each threshold.
    threshvals = (0:(divs-1)) * (ROIMax / (divs - 1));
    
    lr_threshvals = cell(1, numel(threshvals));
    num_threshvals = 0;
    
    for i = 1:numel(threshvals)
        threshval = threshvals(i);
        l_threshvals = lvals_nonnegative(lvals_nonnegative >= threshval);
        r_threshvals = rvals_nonnegative(rvals_nonnegative >= threshval);
        
        %- Check if both hemispheres have enough voxels above the threshold
        % (MIN_NUM_THRESH_VOXELS)
        if (numel(l_threshvals) >= MIN_NUM_THRESH_VOXELS && numel(r_threshvals) >= MIN_NUM_THRESH_VOXELS)
            num_threshvals = num_threshvals + 1;
            lr_threshvals{num_threshvals} = {l_threshvals, r_threshvals};
        else
            break;
        end
    end
    interval_length(j) = length(lr_threshvals);
    
    TB_LIs = zeros(num_threshvals, n_resampling);
    
    for i = 1:num_threshvals
        l_set = lr_threshvals{i}{1};
        r_set = lr_threshvals{i}{2};
        
        l_n = max(round(numel(l_set) * RESAMPLE_RATIO), MIN_NUM_THRESH_VOXELS);
        r_n = max(round(numel(r_set) * RESAMPLE_RATIO), MIN_NUM_THRESH_VOXELS);
        
        l_indices = randi(numel(l_set), n_resampling, l_n);
        lactivity = sum(l_set(l_indices), 2);
        
        r_indices = randi(numel(r_set), n_resampling, r_n);
        ractivity = sum(r_set(r_indices), 2);
        
        TB_LIs(i, :) = (lactivity - ractivity) ./ (lactivity + ractivity);
    end
    
    weights = 1:num_threshvals;
    numerator = sum(trimmean(TB_LIs, 0.25, 2)' .* weights);
    weighted_li(j,:) = (numerator / sum(weights) * 100);
end
% ft_progress('close');

[~, idx_mx] = max(weighted_li); LI_max = wi(idx_mx,:);

if cfg_main.fplot ==1
    
    figure,plot(weighted_li),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1000   300]);
    %     title([cfg_main.tit, ' - ', tmp.Comment]),
    xlabel('temporal windows (sec)')
    ylabel('LI')
    set(gca,'color','none');
end

end
