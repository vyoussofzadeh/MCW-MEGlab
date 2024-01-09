function varargout = process_roi_analysis(varargin )
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
sProcess.Comment     = 'ROI analysis, HCP atlas';
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

% Modify time interval input
sProcess.options.time_interval.Comment = 'Choose a time interval:';
sProcess.options.time_interval.Type    = 'combobox';
sProcess.options.time_interval.Value   = {1, {'Time Interval', 'Averaged Sources'}};

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

% Select effect type (Positive, Negative, Absolute values)
effect = selectEffectType(sProcess.options.effect.Value{1});

% Process ImageGridAmp based on selected effect
ImageGridAmp = processImageGridAmp(sResultP.ImageGridAmp, effect);

%%
ProtocolInfo = bst_get('ProtocolInfo');
SurfaceFile = load(fullfile(ProtocolInfo.SUBJECTS, sResultP.SurfaceFile));

Scouts = [];
sScout = [];
for i = 1:length(SurfaceFile.Atlas)
    if contains(SurfaceFile.Atlas(i).Name, {'mmp_in_mni_symmetrical_1'})
        Scouts = SurfaceFile.Atlas(i).Scouts;
    end
end
sScout.Scouts = Scouts;

% raw value (whole-brain)
rois = []; pow_parcel = [];
for i=1:length(sScout.Scouts)
    pow_parcel(i) = mean(ImageGridAmp(sScout.Scouts(i).Vertices));
end

rois = {sScout.Scouts.Label};
region = {sScout.Scouts.Region};

%%
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';
cfg = struct('src_fname', [], 'glass_dir', glass_dir, 'glass_atlas', sScout, 'plotflag', 0);
Data_hcp_atlas = ecpfunc_hcp_atlas(cfg);

%%
group_members = Data_hcp_atlas.glass_net_L_label;

idx_L = [];
for i=1:length(group_members)
    grois = group_members{i};
    idx = [];
    for j=1:length(grois)
        idx(j) = strmatch(grois{j},rois);
        %         disp([grois{j}, rois(idx(j))])
    end
    idx_L{i} = idx;
end

group_members = Data_hcp_atlas.glass_net_R_label;

idx_R = [];
for i=1:length(group_members)
    grois = group_members{i}(1:end);
    idx = [];
    for j=1:length(grois)
        idx(j) = strmatch(grois{j},rois);
        %         disp([grois{j}, rois(idx(j))])
    end
    idx_R{i} = idx;
end

%%
% Compute power for each group member in the left hemisphere
pow_parcel_L = [];
for i=1:length(idx_L)
    group_roi_indices = idx_L{i};
    group_roi_power = 0;
    for j=1:length(group_roi_indices)
        roi_vertices = Scouts(group_roi_indices(j)).Vertices;
        group_roi_power = group_roi_power + mean(ImageGridAmp(roi_vertices));
    end
    pow_parcel_L(i) = group_roi_power / length(group_roi_indices);
end

% Compute power for each group member in the right hemisphere
pow_parcel_R = [];
for i=1:length(idx_R)
    group_roi_indices = idx_R{i};
    group_roi_power = 0;
    for j=1:length(group_roi_indices)
        roi_vertices = Scouts(group_roi_indices(j)).Vertices;
        group_roi_power = group_roi_power + mean(ImageGridAmp(roi_vertices));
    end
    pow_parcel_R(i) = group_roi_power / length(group_roi_indices);
end

% Optionally, combine the power values for left and right hemispheres
pow_parcel_combined = [pow_parcel_L, pow_parcel_R];

%% Report
clc
% Define group names
group_names = {'Angular', 'Frontal', 'Occipital', 'Other', 'PCingPrecun', ...
    'Temporal', 'BTLA', 'VWFA', 'ATG', 'PSTG', 'lateral'};

% Ensure the number of groups matches the number of names
if length(group_names) ~= length(idx_L) || length(group_names) ~= length(idx_R)
    error('Number of groups in idx_L/idx_R does not match the number of provided group names.');
end

% Create a table for left hemisphere
pow_table_L = table(group_names', pow_parcel_L', 'VariableNames', {'Group', 'Power_Left'});

% Create a table for right hemisphere
pow_table_R = table(group_names', pow_parcel_R', 'VariableNames', {'Group', 'Power_Right'});

LI = 100.*(pow_parcel_L - pow_parcel_R)./(pow_parcel_L + pow_parcel_R);

% Optionally, combine both tables
pow_table_combined = table(group_names', pow_parcel_L', pow_parcel_R', LI', 'VariableNames', {'network_ROIs', 'Power_Left', 'Power_Right', 'LI'});

disp(pow_table_combined);



%%
% Assuming rois and pow_parcel are already defined
% rois - a cell array containing the labels of each ROI
% pow_parcel - an array containing the power values for each ROI

% Check if the number of ROI labels matches the number of power values
if length(rois) ~= length(pow_parcel)
    error('The number of ROI labels does not match the number of power values.');
end

roi_table = table(rois', pow_parcel', region', 'VariableNames', {'ROI_Label', 'Power', 'Region'});

%%
% Assuming roi_table is already defined with columns 'ROI_Label' and 'Power'

% Number of top entries to display
numTopEntries = 10;

% Sort the table by 'Power' in descending order
sorted_table = sortrows(roi_table, 'Power', 'descend');

% Select the top entries
top_entries = sorted_table(1:min(end, numTopEntries), :);

% Display the top entries
disp(top_entries);

%%
% Assuming rois, pow_parcel, and region are already defined
% rois - a cell array containing the labels of each ROI
% pow_parcel - an array containing the power values for each ROI
% region - a cell array containing the region for each ROI

% Check if the number of elements match
if length(rois) ~= length(pow_parcel) || length(rois) ~= length(region)
    error('The number of elements in rois, pow_parcel, and region must be the same.');
end

% Create a table with ROI labels, their corresponding power values, and regions
roi_table = table(rois', pow_parcel', region', 'VariableNames', {'ROI_Label', 'Power', 'Region'});

% Get unique regions
unique_regions = unique(region);

% Create a table to store the aggregated results
aggregated_table = table('Size', [0 3], 'VariableTypes', {'string', 'double', 'string'}, 'VariableNames', {'Region', 'Average_Power', 'ROIs'});

% Iterate over each unique region and calculate the average power
for i = 1:length(unique_regions)
    current_region = unique_regions{i};
    current_group_indices = strcmp(roi_table.Region, current_region);
    average_power = mean(roi_table.Power(current_group_indices));
    rois_in_group = roi_table.ROI_Label(current_group_indices);
    aggregated_table = [aggregated_table; {current_region, average_power, strjoin(rois_in_group, ', ')}];
end

% Display the aggregated table
disp(aggregated_table);

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



function Data_hcp_atlas = ecpfunc_hcp_atlas(cfg_main)
% ECP functions
% Project: ECP_SD
% Written by: Vahab Youssof Zadeh
% Update: 05/31/2023

src_fname = cfg_main.src_fname;
glass_dir = cfg_main.glass_dir;
glass_atlas = cfg_main.glass_atlas;

%%
% '/group/jbinder/ECP/MEG/laterality_index/bilateral_glasser_lateral.tsv'

%%
% Load atlas
% atlas = load('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat');
% atlas = load(glass_atlas);
atlas = (glass_atlas);
groups_labels = {'Angular', 'Frontal', 'Occipital', 'Other', 'PCingPrecun', 'Temporal'};

rois = {atlas.Scouts.Label};
region = {atlas.Scouts.Region};

all_idx_L = find(startsWith(rois, 'L_'))';
all_idx_R = find(startsWith(rois, 'R_'))';

if cfg_main.plotflag == 1
    
    % Plot HCP atlas
    cfg = struct();
    cfg.atlas = atlas;
    cfg.src_fname = src_fname;
    % '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
    cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
    cfg.index_L = all_idx_L;
    cfg.index_R = all_idx_R;
    cfg.rois = rois;
    cfg.rois_sel = 1:180;
    cfg.title = '';
    do_plot_HCP_atlas(cfg)
    
end

% Update frontal region
idx_sel_L = strcmp(region(all_idx_L), 'LF');
idx_sel_R = strcmp(region(all_idx_R), 'RF');

% Load atlas ROIs from fMRI study
% '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/Glasser';
load(fullfile(glass_dir, 'LI_glasser_manual_net_12.mat'), 'glass_net_L_label', 'glass_net_R_label');

% Update frontal region labels
glass_net_L_label{2} = [glass_net_L_label{2}; rois(all_idx_L(idx_sel_L))'];
glass_net_R_label{2} = [glass_net_R_label{2}; rois(all_idx_R(idx_sel_R))'];

% Add BTLA labels
btla = [2, 3, 5, 8, 9, 16, 17, 18, 21, 22]; net_sel = 6;
% glass_net_L_label{6} = glass_net_L_label{6}(btla);
% glass_net_R_label{6} = glass_net_R_label{6}(btla);

BTLA_L_label = [];
BTLA_R_label = [];
for i=1:length(btla)
    BTLA_L_label = [BTLA_L_label; glass_net_L_label{net_sel}(btla(i))];
    BTLA_R_label = [BTLA_R_label; glass_net_R_label{net_sel}(btla(i))];
end

glass_net_L_label{7} = BTLA_L_label;
glass_net_R_label{7} = BTLA_R_label;

groups_labels{7} = 'BTLA';

%% Add VWFA labels
vw2 = [6, 14, 15, 81]; net_sel = 4;
% glass_net_L_label{4} = glass_net_L_label{4}(vw2);
% glass_net_R_label{4} = glass_net_R_label{4}(vw2);

VW_L_label = [];
VW_R_label = [];
for i=1:length(vw2)
    VW_L_label = [VW_L_label; glass_net_L_label{net_sel}(vw2(i))];
    VW_R_label = [VW_R_label; glass_net_R_label{net_sel}(vw2(i))];
end

glass_net_L_label{8} = VW_L_label;
glass_net_R_label{8} = VW_R_label;

groups_labels{8} = 'VWFA';

% Add LT and RT region labels
% idx_sel_L = strcmp(region(all_idx_L), 'LT');
% idx_sel_R = strcmp(region(all_idx_R), 'RT');
% glass_net_L_label{7} = rois(all_idx_L(idx_sel_L));
% glass_net_R_label{7} = rois(all_idx_R(idx_sel_R));


%% Add ATG labels
ATG_labels = {'L_TGv_ROI', 'L_TGd_ROI'};

% Find indices that correspond to ATG and STG labels
ATG_L_label = [];
ATG_R_label = [];
for i=1:length(ATG_labels)
    ATG_L_label = [ATG_L_label; rois(strcmp(rois, ATG_labels{i}))];
    ATG_R_label = [ATG_R_label; rois(strcmp(rois, strrep(ATG_labels{i}, 'L_', 'R_')))];
end

glass_net_L_label{9} = ATG_L_label;
glass_net_R_label{9} = ATG_R_label;

groups_labels{9} = 'ATG'; %Anterior temporal G.

%% Post. STG
PSTG_labels = {'L_TE1p_ROI'};

% Find indices that correspond to ATG and STG labels
PSTG_L_label = [];
PSTG_R_label = [];
for i=1:length(PSTG_labels)
    PSTG_L_label = [PSTG_L_label; rois(strcmp(rois, PSTG_labels{i}))];
    PSTG_R_label = [PSTG_R_label; rois(strcmp(rois, strrep(PSTG_labels{i}, 'L_', 'R_')))];
end

glass_net_L_label{10} = PSTG_L_label;
glass_net_R_label{10} = PSTG_R_label;

groups_labels{10} = 'PSTG'; %Post. STG


%% Lateral rois
% see, /data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Pipe_check_atlas.m for details
LI_glasser_lateral_rois = load(fullfile(glass_dir,'LI_glasser_lateral_rois.mat'));
groups_labels = [groups_labels, 'lateral'];

glass_net_L_label{11} = LI_glasser_lateral_rois.glass_roi_lat_L_name;
glass_net_R_label{11} = LI_glasser_lateral_rois.glass_roi_lat_R_name;


%%
Data_hcp_atlas.glass_net_L_label = glass_net_L_label;
Data_hcp_atlas.glass_net_R_label = glass_net_R_label;
Data_hcp_atlas.groups_labels = groups_labels;
Data_hcp_atlas.atlas = atlas;
Data_hcp_atlas.rois = rois;

end