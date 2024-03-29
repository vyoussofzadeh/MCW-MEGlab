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

% Add options to select the LI computation method as checkboxes
sProcess.options.methodCounting.Comment = 'Use Counting Method';
sProcess.options.methodCounting.Type    = 'checkbox';
sProcess.options.methodCounting.Value   = 0; % Not selected by default

sProcess.options.methodBootstrap.Comment = 'Use Bootstrapping Method';
sProcess.options.methodBootstrap.Type    = 'checkbox';
sProcess.options.methodBootstrap.Value   = 0; % Not selected by default

% Add bootstrap parameters
sProcess.options.divs.Comment = 'Number of divisions for bootstrap:';
sProcess.options.divs.Type    = 'value';
sProcess.options.divs.Value   = {10, '', 1, 100, 1}; % default 10, min 1, max 100, step 1

sProcess.options.n_resampling.Comment = 'Number of resampling iterations:';
sProcess.options.n_resampling.Type    = 'value';
sProcess.options.n_resampling.Value   = {20, '', 1, 1000, 1}; % default 20, min 1, max 1000, step 1

sProcess.options.RESAMPLE_RATIO.Comment = 'Resample ratio (%):';
sProcess.options.RESAMPLE_RATIO.Type    = 'value';
sProcess.options.RESAMPLE_RATIO.Value   = {75, '%', 0, 100, 1, 1};  % Now represented as a percentage

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

% Handle the selected LI computation method
if sProcess.options.methodCounting.Value ==1
    computeLI(cfg_LI);  % Counting-based method
end
if sProcess.options.methodBootstrap.Value ==1
    cfg_LI.divs =  sProcess.options.divs.Value{1};  % Adjust as needed
    cfg_LI.n_resampling =  sProcess.options.n_resampling.Value{1};  % Adjust as needed
    cfg_LI.RESAMPLE_RATIO = sProcess.options.RESAMPLE_RATIO.Value{1} / 100;  % Adjust as needed
    computeLI_bootstrap(cfg_LI);
end

if time_interval ~=2
    disp(['Time: ', num2str(timerange(1)), '-', num2str(timerange(2))]);
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

expectedRegions = {...
'bankssts L', 'bankssts R', ...
'caudalanteriorcingulate L', 'caudalanteriorcingulate R', ...
'caudalmiddlefrontal L', 'caudalmiddlefrontal R', ...
'cuneus L', 'cuneus R', ...
'entorhinal L', 'entorhinal R', ...
'frontalpole L', 'frontalpole R', ...
'fusiform L', 'fusiform R', ...
'inferiorparietal L', 'inferiorparietal R', ...
'inferiortemporal L', 'inferiortemporal R', ...
'insula L', 'insula R', ...
'isthmuscingulate L', 'isthmuscingulate R', ...
'lateraloccipital L', 'lateraloccipital R', ...
'lateralorbitofrontal L', 'lateralorbitofrontal R', ...
'lingual L', 'lingual R', ...
'medialorbitofrontal L', 'medialorbitofrontal R', ...
'middletemporal L', 'middletemporal R', ...
'paracentral L', 'paracentral R', ...
'parahippocampal L', 'parahippocampal R', ...
'parsopercularis L', 'parsopercularis R', ...
'parsorbitalis L', 'parsorbitalis R', ...
'parstriangularis L', 'parstriangularis R', ...
'pericalcarine L', 'pericalcarine R', ...
'postcentral L', 'postcentral R', ...
'posteriorcingulate L', 'posteriorcingulate R', ...
'precentral L', 'precentral R', ...
'precuneus L', 'precuneus R', ...
'rostralanteriorcingulate L', 'rostralanteriorcingulate R', ...
'rostralmiddlefrontal L', 'rostralmiddlefrontal R', ...
'superiorfrontal L', 'superiorfrontal R', ...
'superiorparietal L', 'superiorparietal R', ...
'superiortemporal L', 'superiortemporal R', ...
'supramarginal L', 'supramarginal R', ...
'temporalpole L', 'temporalpole R', ...
'transversetemporal L', 'transversetemporal R'...
};

% Handle case when number of anatomical regions are not identical to atlas regions
actualRegions = {sScout.Scouts.Label};
missingRegions = setdiff(expectedRegions, actualRegions);

% Assuming sScout.Scouts is not empty and has at least one scout
if ~isempty(sScout.Scouts)
    % Identify all fields from the first scout as a template
    fieldNames = fieldnames(sScout.Scouts(1));
    % Prepare an empty scout template with all fields
    emptyScout = cell2struct(cell(length(fieldNames), 1), fieldNames, 1);
    
    % Default empty values for known fields
    emptyScout.Label = ''; % Update as necessary
    emptyScout.Vertices = [];
    emptyScout.Seed = 0; % Or any appropriate 'empty' value
    
    % For any other fields in your structure, set an appropriate 'empty' value
    % For example, if there's a 'Color' field, you might do:
    % emptyScout.Color = [0, 0, 0]; % Assuming color is RGB
    % Adjust the above line for each additional field with a sensible 'empty' value
    
    % Now, handle missing regions with this updated emptyScout
    for i = 1:length(missingRegions)
        emptyScout.Label = missingRegions{i};
        % Insert the empty scout at the correct position
        idx = find(strcmp(expectedRegions, missingRegions{i}));
        sScout.Scouts = [sScout.Scouts(1:idx-1), emptyScout, sScout.Scouts(idx:end)];
    end
else
    warning('sScout.Scouts is empty, cannot determine structure fields.');
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
% Computes the Laterality Index (LI) using bootstrapping and exports the results, including vertex counts.

% Perform bootstrapping for each ROI
RoiLabels = cfg_LI.RoiLabels;
TotROI = length(cfg_LI.RoiIndices);
Summ_LI = zeros(1, TotROI); % Initialize the vector for final LIs
L_vertices_total = zeros(1, TotROI); % Initialize vectors for left vertices count
R_vertices_total = zeros(1, TotROI); % Initialize vectors for right vertices count
LI_label_out = cell(1, TotROI);

disp('Bootstrapping ..')
for ii = 1:TotROI
    
    disp(['Processing ROI ', num2str(ii), ' of ', num2str(TotROI)])
    cfg_main = [];
    cfg_main.atlas = cfg_LI.sScout;
    cfg_main.RoiIndices = cfg_LI.RoiIndices{ii}; % Pass current ROI indices
    cfg_main.divs = cfg_LI.divs; % Adjust as needed
    cfg_main.n_resampling = cfg_LI.n_resampling; % Corrected to use n_resampling from cfg_LI
    cfg_main.RESAMPLE_RATIO = cfg_LI.RESAMPLE_RATIO; % Adjust as needed
    cfg_main.t1 = cfg_LI.t1;
    cfg_main.t2 = cfg_LI.t2;
    cfg_main.ImageGridAmp = cfg_LI.ImageGridAmp;
    cfg_main.time_interval = cfg_LI.time_interval;
    
    % Call bootstrapping function for current ROI
    [weighted_li, ~, L_vertices_above_thresh, R_vertices_above_thresh] = do_LI_bootstrap(cfg_main);
    
    % Store results
    Summ_LI(ii) = weighted_li; % Assuming weighted_li represents the LI for simplicity
    L_vertices_total(ii) = sum(L_vertices_above_thresh); % Summing vertices counts for left hemisphere
    R_vertices_total(ii) = sum(R_vertices_above_thresh); % Summing vertices counts for right hemisphere
    LI_label_out{ii} = RoiLabels{ii};
end

% Save or display results
savedir = cfg_LI.savedir; % Directory to save the results
sname = cfg_LI.sname; % Use sname from cfg_LI for file naming
filename = fullfile(savedir, sname);

% Open file for writing
fid = fopen(filename, 'w');
fprintf(fid, 'ROI\tLI\tL_Vertices\tR_Vertices\n'); % Updated header to include vertex counts
for i = 1:TotROI
    fprintf(fid, '%s\t%f\t%d\t%d\n', LI_label_out{i}, Summ_LI(i), L_vertices_total(i), R_vertices_total(i));
end
fclose(fid);

disp(['Results saved to: ', filename]);

% Optional: Display results in console for quick viewing
disp('Bootstrap LI Results with Vertex Counts:');
disp(table(RoiLabels', Summ_LI', L_vertices_total', R_vertices_total', 'VariableNames', {'ROI', 'LI', 'L_Vertices', 'R_Vertices'}));

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

function [weighted_li, num_threshvals, L_vertices_above_thresh, R_vertices_above_thresh] = do_LI_bootstrap(cfg_main)
% Modified to use RoiIndices for distinguishing between left and right hemisphere ROIs
% and to refine the reporting of vertices above threshold across bootstrapping.

% Extract necessary configurations
divs = cfg_main.divs;
n_resampling = cfg_main.n_resampling;
RESAMPLE_RATIO = cfg_main.RESAMPLE_RATIO;
RoiIndices = cfg_main.RoiIndices;
MIN_NUM_THRESH_VOXELS = divs / RESAMPLE_RATIO; % Adjust based on your specific logic


if cfg_main.time_interval == 2 || size(cfg_main.ImageGridAmp,2) == 1
    ImageGridAmp = mean(cfg_main.ImageGridAmp,2);
else
    ImageGridAmp = cfg_main.ImageGridAmp(:, cfg_main.t1:cfg_main.t2);
end

if size(ImageGridAmp,2) > 500   
    ImageGridAmp = ImageGridAmp(:,1:10:end);    
end

% Initialize output variables
weighted_li = 0;
num_threshvals = 0;

% Initialize cumulative variables for counting vertices above threshold
cumulative_L_vertices_above_thresh = 0;
cumulative_R_vertices_above_thresh = 0;

% Distinguish between left and right hemisphere ROIs
LHscout = horzcat(cfg_main.atlas.Scouts(RoiIndices(1:2:end)).Vertices);
RHscout = horzcat(cfg_main.atlas.Scouts(RoiIndices(2:2:end)).Vertices);

% Extract amplitude values for left and right regions
LHvals = ImageGridAmp(LHscout, :);
RHvals = ImageGridAmp(RHscout, :);

% Determine thresholds based on max value across both hemispheres
ROIMax = max([LHvals(:); RHvals(:)]);
threshvals = linspace(0, ROIMax, divs);

if divs ==1
    disp('divs should be greater than 1.')
end

% Perform bootstrapping for each threshold
for thresh_idx = 1:numel(threshvals)
    
    thresh = threshvals(thresh_idx);
    l_above_thresh = LHvals >= thresh;
    r_above_thresh = RHvals >= thresh;
    
    % Ensure both hemispheres have enough voxels above threshold
    if sum(l_above_thresh(:)) < MIN_NUM_THRESH_VOXELS && sum(r_above_thresh(:)) < MIN_NUM_THRESH_VOXELS
        break; % Stop if not enough voxels
    end
    
    cumulative_L_vertices_above_thresh = cumulative_L_vertices_above_thresh + sum(l_above_thresh(:));
    cumulative_R_vertices_above_thresh = cumulative_R_vertices_above_thresh + sum(r_above_thresh(:));
    
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

% Adjust final vertex counts to reflect average contribution if needed
L_vertices_above_thresh = round(cumulative_L_vertices_above_thresh / num_threshvals);
R_vertices_above_thresh = round(cumulative_R_vertices_above_thresh / num_threshvals);

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
