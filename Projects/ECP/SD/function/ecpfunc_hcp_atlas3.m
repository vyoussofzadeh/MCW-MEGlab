function Data_hcp_atlas = ecpfunc_hcp_atlas3(cfg_main)

% ECP functions
% Project: ECP_SD
% Written by: Vahab Youssof Zadeh
% Update: 05/31/2023

src_fname = cfg_main.src_fname;
glass_dir = cfg_main.glass_dir;
glass_atlas = cfg_main.glass_atlas;

%%
atlas = load(glass_atlas);
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
load(fullfile(glass_dir, 'LI_glasser_manual_net_12.mat'), 'glass_net_L_label', 'glass_net_R_label');

% Update frontal region labels
glass_net_L_label{2} = unique([glass_net_L_label{2}; rois(all_idx_L(idx_sel_L))']);
glass_net_R_label{2} = unique([glass_net_R_label{2}; rois(all_idx_R(idx_sel_R))']);

%% Add BTLA labels
btla = [2, 3, 5, 8, 9, 16, 17, 18, 21, 22]; net_sel = 6;

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

VW_L_label = [];
VW_R_label = [];
for i=1:length(vw2)
    VW_L_label = [VW_L_label; glass_net_L_label{net_sel}(vw2(i))];
    VW_R_label = [VW_R_label; glass_net_R_label{net_sel}(vw2(i))];
end

glass_net_L_label{8} = VW_L_label;
glass_net_R_label{8} = VW_R_label;

groups_labels{8} = 'VWFA';

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
% see,
% /data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Pipe_check_atlas.m
% for details - OLD

% filePath = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser/bilateral_glasser_lateral.tsv';
%
% % Specify the file path (adjust this to your local setup)
% % filePath = 'path_to_your_file/bilateral_glasser_lateral.tsv';
%
% % Read the TSV file into a table
% data = readtable(filePath, 'FileType', 'text', 'Delimiter', '\t');
%
% % Filter rows where roi_num is 2 (indicating lateral)
% lateralROIs = data(data.roi_num == 2, :);
%
% % Generate corresponding left parcel and roi names by replacing 'R_' with 'L_'
% lateralROIs.L_parcel_name = replace(lateralROIs.parcel_name, 'R_', 'L_');
% lateralROIs.L_roi_name = replace(lateralROIs.roi_name, 'R_', 'L_');
%
% % Display the relevant columns
% disp(lateralROIs(:, {'parcel_name', 'roi_name', 'L_parcel_name', 'L_roi_name'}));

LI_glasser_lateral_rois = load(fullfile(glass_dir,'LI_glasser_lateral_rois.mat'));
groups_labels = [groups_labels, 'lateral'];

glass_net_L_label{11} = LI_glasser_lateral_rois.glass_roi_lat_L_name;
glass_net_R_label{11} = LI_glasser_lateral_rois.glass_roi_lat_R_name;

%% Initialize containers for network ROI indices
network_roi_indices_L = cell(1, length(groups_labels));
network_roi_indices_R = cell(1, length(groups_labels));

%% Iterate through each network and find/store indices
for i = 1:length(groups_labels)
    network_roi_indices_L{i} = find(ismember(rois(all_idx_L), glass_net_L_label{i}));
    network_roi_indices_R{i} = find(ismember(rois(all_idx_R), glass_net_R_label{i}));
end

%%
Data_hcp_atlas.glass_net_L_label = glass_net_L_label;
Data_hcp_atlas.glass_net_R_label = glass_net_R_label;
Data_hcp_atlas.groups_labels = groups_labels;
Data_hcp_atlas.atlas = atlas;
Data_hcp_atlas.rois = rois;
Data_hcp_atlas.network_roi_indices_L = network_roi_indices_L;
Data_hcp_atlas.network_roi_indices_R = network_roi_indices_R;

end