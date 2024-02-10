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
sProcess.Comment     = 'Compute LI, surface-based, DK atlas';
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

% Modify time interval input to include "Averaged Time Interval"
sProcess.options.time_interval.Comment = 'Choose a time interval:';
sProcess.options.time_interval.Type    = 'combobox';
sProcess.options.time_interval.Value   = {1, {'Specific Time Interval', 'Averaged Time Interval'}};

% Active time window
sProcess.options.poststim.Comment = 'Enter specific time interval:';
sProcess.options.poststim.Type    = 'poststim';
sProcess.options.poststim.Value   = [];

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

% Conditional activation of specific time interval input
if time_interval == 1
    timerange = sProcess.options.poststim.Value{1};
else
    timerange =1; % Follow the existing procedure for other options
end

% Select effect type (Positive, Negative, Absolute values)
effect = selectEffectType(sProcess.options.effect.Value{1});

% Define ROI-related parameters

Ratio4Threshold = sProcess.options.ratio4threshold.Value{1}/100;

% Define threshold type
Threshtype = selectThresholdType(sProcess.options.threshtype.Value{1});

% Process ImageGridAmp based on selected effect
ImageGridAmp = processImageGridAmp(sResultP.ImageGridAmp, effect);

% Determine max values for various time windows
[AllMax, GlobalMax, t1, t2] = determineMaxValues(time_interval, ImageGridAmp, sResultP, timerange);

% Convert Desikan-Killiany scout to selected scouts
[sScout, ~] = convertDesikanKillianyScout(sResultP);

% Example modification in the part that processes ImageGridAmp
if sProcess.options.time_interval.Value{1} == 2 && length(sResultP.Time) > 10
    timerange = sProcess.options.poststim.Value{1};
    % Compute the average across the selected time range
    [~, t1] = min(abs(sResultP.Time - timerange(1)));
    [~, t2] = min(abs(sResultP.Time - timerange(2)));
    ImageGridAmp = mean(ImageGridAmp(:, t1:t2), 2);
else
    % Existing logic for specific time intervals or other options
end

% Define ROIs
[RoiLabels, RoiIndices] = defineROIs();

% Compute LI
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
        computeLI(cfg_LI);
    case 2
        computeLI_bootstrap(cfg_LI)
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

function effect = selectEffectType(effect)
% Prompt user to select effect type

% Ensure valid selection
if isempty(effect) || ~any(effect == [1, 2, 3])
    error('Invalid effect type selection. Choose 1, 2, or 3.');
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
    [~, IA] = setdiff(rois_temp, rois);
    
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

function computeLI_bootstrap(cfg_LI)
% Computes the Laterality Index (LI) using bootstrapping and exports the results.

% Setup - Assuming cfg_LI contains all necessary configurations including:
% ImageGridAmp, Time, RoiLabels, RoiIndices, sScout, t1, t2, savedir

% Perform bootstrapping for each ROI
RoiLabels = cfg_LI.RoiLabels;
TotROI = length(cfg_LI.RoiIndices);
Summ_LI = zeros(1, TotROI); % Initialize the vector for final LIs
L_count = zeros(1, TotROI);
R_count = zeros(1, TotROI);
LI_label_out = cell(1, TotROI);

for ii = 1:TotROI
    
    disp(['network ', num2str(ii), ' of ', num2str(TotROI)])
    cfg_main = [];
    cfg_main.atlas = cfg_LI.sScout;
    cfg_main.RoiIndices = cfg_LI.RoiIndices{ii}; % Pass current ROI indices
    cfg_main.divs = 10; % Example, adjust as needed
    cfg_main.n_resampling = 100; % Example, adjust as needed
    cfg_main.RESAMPLE_RATIO = 0.75; % Example, adjust as needed
    cfg_main.t1 = cfg_LI.t1;
    cfg_main.t2 = cfg_LI.t2;
    cfg_main.ImageGridAmp = cfg_LI.ImageGridAmp;
    
    % Call bootstrapping function for current ROI
    [weighted_li, ~] = do_LI_bootstrap(cfg_main);
    
    % Store results
    Summ_LI(ii) = weighted_li; % Assuming weighted_li represents the LI for simplicity
    LI_label_out{ii} = RoiLabels{ii};
    
end

% Save or display results
savedir = cfg_LI.savedir; % Directory to save the results
sname = 'bootstrap_results.xls'; % Example filename, adjust as needed
filename = fullfile(savedir, sname);

% Open file for writing
fid = fopen(filename, 'w');
fprintf(fid, 'ROI\tLI\tLeft_count\tRight_count\n');
for i = 1:TotROI
    fprintf(fid, '%s\t%f\t%d\t%d\n', LI_label_out{i}, Summ_LI(i), L_count(i), R_count(i));
end
fclose(fid);

disp(['Results saved to: ' filename]);

%%
% Display the path to the saved file
disp(['Results saved to: ' filename]);

% Optional: Display results
disp('=================')
disp('                 ')
a = table(RoiLabels'); a.Properties.VariableNames{'Var1'} = 'ROI';
b = table(Summ_LI'); b.Properties.VariableNames{'Var1'} = 'LI';
disp([a,b])

end

function computeLI(cfg_LI)
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

TotROI = length(RoiIndices);

s1='LI_';
Summ_LI=zeros(1,TotROI); % initialize the vector that summarizes the final LIs  % added JL 11212014
Summ_LI_Label='ROI Labels: '; % initialize the string that summarizes the ROI labels  % added JL 11212014

% Adjustments for dimensions of verticies when concatenating
for i=1:length(sScout.Scouts)
    % Check if Vertices is nx1, transpose only in this case
    if size(sScout.Scouts(i).Vertices, 1) > 1
        sScout.Scouts(i).Vertices = sScout.Scouts(i).Vertices';
    end
end

%%
LI_label_out={};
for ii = 1:length(RoiIndices)
    
    s2 = RoiLabels{ii};
    hemi_roi_num=length(RoiIndices{ii});
    curr_subregion=sScout.Scouts(RoiIndices{ii});
    
    %Odd indices are left Rois
    Ltemp_region = [];
    Ltemp_label  = [];
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
    L_ROIcount = length(ind_L);
    L_count(ii) = L_ROIcount;
    %JL@10/30/14    indHalfROImaxR = find(RHvals > ROIMax/2); % Not being used right now
    R_ROIcount = length(ind_R);
    R_count(ii) = R_ROIcount;
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

%%
tempfile = fopen(filename, 'w');

% Print headers
fprintf(tempfile, 'ROI\tLI\tLeft_count\tRight_count\n');

% Ensure L_count and R_count are column vectors
L_count1 = L_count(:);
R_count1 = R_count(:);

% Printing the labels, LI values, L_count, and R_count on separate lines
for i = 1:length(LI_label_out)
    fprintf(tempfile, '%s\t', LI_label_out{i});
    fprintf(tempfile, '%f\t', Summ_LI(i));
    fprintf(tempfile, '%d\t', L_count1(i));
    fprintf(tempfile, '%d\n', R_count1(i));
end

% Adding an extra newline for separation
fprintf(tempfile, '\n');

% Printing the threshold
fprintf(tempfile, 'Threshold\t%f\n', threshold);

fclose(tempfile);

%%
% Display the path to the saved file
disp(['Results saved to: ' filename]);

% Optional: Display results
disp('=================')
disp('                 ')
a = table(RoiLabels'); a.Properties.VariableNames{'Var1'} = 'ROI';
b = table(Summ_LI'); b.Properties.VariableNames{'Var1'} = 'LI';
c = table([L_count;R_count]'); c.Properties.VariableNames{'Var1'} = 'Left_vs_right';
d = [a,b,c];
disp(d)

end

function [weighted_li, num_threshvals] = do_LI_bootstrap(cfg_main)
% Modified to use RoiIndices for distinguishing between left and right hemisphere ROIs.

% Extract necessary configurations
divs = cfg_main.divs;
n_resampling = cfg_main.n_resampling;
RESAMPLE_RATIO = cfg_main.RESAMPLE_RATIO;
% time_interval = [cfg_main.t1, cfg_main.t2];
RoiIndices = cfg_main.RoiIndices; % Assuming RoiIndices is passed in cfg_main
MIN_NUM_THRESH_VOXELS = round(5 / RESAMPLE_RATIO);

% Assuming the ImageGridAmp data is already in cfg_main
ImageGridAmp = cfg_main.ImageGridAmp;

% Initialize output variables
weighted_li = 0;
num_threshvals = 0;

% Preprocess data for the specified time interval
% d_in = mean(ImageGridAmp(:, time_interval(1):time_interval(2)), 2);
% d_in = ImageGridAmp;
% ImageGridAmp = (d_in); % Ensure non-negative values

%%
hemi_roi_num=length(RoiIndices);
curr_subregion=cfg_main.atlas.Scouts(RoiIndices);

%Odd indices are left Rois
Ltemp_region = [];
Ltemp_label  = [];
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

%%

% Extract amplitude values for left and right regions
LHvals = ImageGridAmp(LHscout,:);
RHvals = ImageGridAmp(RHscout,:);

% Determine thresholds based on max value across both hemispheres
ROIMax = max([max(LHvals), max(RHvals)]);
threshvals = linspace(0, ROIMax, divs);

% Perform bootstrapping for each threshold
for thresh_idx = 1:numel(threshvals)
    thresh = threshvals(thresh_idx);
    l_above_thresh = LHvals(LHvals >= thresh);
    r_above_thresh = RHvals(RHvals >= thresh);
    
    % Ensure both hemispheres have enough voxels above threshold
    if numel(l_above_thresh) < MIN_NUM_THRESH_VOXELS || numel(r_above_thresh) < MIN_NUM_THRESH_VOXELS
        break; % Stop if not enough voxels
    end
    
    % Resample and compute LI for the current threshold
    TB_LIs = bootstrapLI(l_above_thresh, r_above_thresh, n_resampling, RESAMPLE_RATIO);
    
    % Weight the LI by threshold index (heavier weight for higher thresholds)
    weighted_li = weighted_li + mean(TB_LIs) * thresh_idx;
    num_threshvals = thresh_idx;
end

% Normalize the weighted LI by the sum of threshold indices
if num_threshvals > 0
    weighted_li = (weighted_li / sum(1:num_threshvals)) * 100;
end

end

function TB_LIs = bootstrapLI(Lvals, Rvals, n_samples, resample_ratio)
% Bootstrap Laterality Index computation for given left and right hemisphere values.
l_n = max(round(numel(Lvals) * resample_ratio), 1);
r_n = max(round(numel(Rvals) * resample_ratio), 1);

TB_LIs = zeros(n_samples, 1);
for i = 1:n_samples
    l_sample = datasample(Lvals, l_n);
    r_sample = datasample(Rvals, r_n);
    TB_LIs(i) = (mean(l_sample) - mean(r_sample)) / (mean(l_sample) + mean(r_sample));
end
end