function Data_hcp_atlas = ecpfunc_hcp_atlas(~)
% ECP functions
% Project: ECP_SD
% Written by: Vahab Youssof Zadeh
% Update: 05/31/2023

% Load atlas
atlas = load('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat');
groups_labels = {'Angular', 'Frontal', 'Occipital', 'Other', 'PCingPrecun', 'Temporal'};

rois = {atlas.Scouts.Label};
region = {atlas.Scouts.Region};

all_idx_L = find(startsWith(rois, 'L_'))';
all_idx_R = find(startsWith(rois, 'R_'))';

% Plot HCP atlas
cfg = struct();
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.index_L = all_idx_L;
cfg.index_R = all_idx_R;
cfg.rois = rois;
cfg.rois_sel = 1:180;
cfg.title = '';
do_plot_HCP_atlas(cfg)

% Update frontal region
idx_sel_L = strcmp(region(all_idx_L), 'LF');
idx_sel_R = strcmp(region(all_idx_R), 'RF');

% Load atlas ROIs from fMRI study
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/Glasser';
load(fullfile(glass_dir, 'LI_glasser_manual_net_12.mat'), 'glass_net_L_label', 'glass_net_R_label');

% Update frontal region labels
glass_net_L_label{2} = [glass_net_L_label{2}; rois(all_idx_L(idx_sel_L))'];
glass_net_R_label{2} = [glass_net_R_label{2}; rois(all_idx_R(idx_sel_R))'];

% Update BTLA labels
btla = [2, 3, 5, 8, 9, 16, 17, 18, 21, 22];
glass_net_L_label{6} = glass_net_L_label{6}(btla);
glass_net_R_label{6} = glass_net_R_label{6}(btla);

% Update VWFA labels
vw2 = [6, 14, 15, 81];
glass_net_L_label{4} = glass_net_L_label{4}(vw2);
glass_net_R_label{4} = glass_net_R_label{4}(vw2);

% Update LT and RT region labels
idx_sel_L = strcmp(region(all_idx_L), 'LT');
idx_sel_R = strcmp(region(all_idx_R), 'RT');
glass_net_L_label{7} = rois(all_idx_L(idx_sel_L));
glass_net_R_label{7} = rois(all_idx_R(idx_sel_R));

Data_hcp_atlas.glass_net_L_label = glass_net_L_label;
Data_hcp_atlas.glass_net_R_label = glass_net_R_label;
Data_hcp_atlas.groups_labels = groups_labels;
Data_hcp_atlas.atlas = atlas;
Data_hcp_atlas.rois = rois;

end
