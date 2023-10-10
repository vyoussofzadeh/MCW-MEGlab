%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 06/27/2023

clear; clc, close('all'); warning off,

%% analysis flags
flags = [];
flags.sourcemaps_contrast = 0;

%% Paths
restoredefaultpath
addpath('./run')
addpath('./function')
Run_setpath

%% Run Brainstorm
Run_BS

%%
Run_load_surface_template

%% THIS PART should be run one time
if flags.sourcemaps_contrast ==1
    cfg = [];
    cfg.protocol = protocol;
    cfg.datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';
    cfg.BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_all_subjects';
    ecpfunc_run_sourcemaps_contrast;
end

%%
clc
cfg = []; 
cfg.protocol = protocol;
cfg.datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';
cfg.BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_all_subjects';
S_data = ecpfunc_read_sourcemaps_contrast(cfg);

%% Subject demog details
clc
cfg = []; cfg.subjs = S_data.subjs; 
cfg.sFiles = S_data.sFiles_32; 
sub_demog_data = ecpfunc_read_sub_demog_contrast(cfg);

%% TLE side (PT only)
patn_neuropsych_data = ecpfunc_read_patn_neuropsych();

%%
clc
cfg = [];
cfg.sub_demog_data = sub_demog_data;
cfg.patn_neuropsych_data = patn_neuropsych_data;
sub_TLE_sub_data = ecpfunc_read_sub_TLE_sub_contrast(cfg);

%%
src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/';
glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat';
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';

%% HCP Atlas
clc, close all
cfg = []; 
cfg.src_fname = src_fname;
cfg.glass_dir = glass_dir;
cfg.glass_atlas = glass_atlas;
Data_hcp_atlas = ecpfunc_hcp_atlas2(cfg);

%% Time intervals (window)
cfg.strt = 0;
cfg.spt = 2;
cfg.overlap = 0.01;
cfg.linterval = 0.3;
wi  = do_time_intervals(cfg);

%%
thre = 0.5;

%%
cfg = [];
cfg.src_fname = src_fname;
cfg.network_sel = [1,2,6];
cfg.Data_hcp_atlas = Data_hcp_atlas;
[idx_L, idx_R, src]  = do_plot_hcp_network(cfg);
net_rois = 'ftp';

%%
disp(Data_hcp_atlas.groups_labels)

cfg = []; cfg.idx_L = idx_L; cfg.idx_R = idx_R; cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.src_fname = src_fname;

cfg.export = 1; cfg.savedir = fullfile(outdir,'group');
cfg.network_sel = [1,2,6]; do_map_HCP_net_sel(cfg);
net_label = 'Fronto_tempro_pri';

% Inspecting atlas areas.
for i=1:11
    cfg.network_sel = i; do_map_HCP_net_sel(cfg);title(Data_hcp_atlas.groups_labels{cfg.network_sel})
end

%%
disp('1: threshold')
disp('2: counting')
disp('3: bootstrapping')
LI_method = input('LI_method sel:');
switch LI_method
    case 1
        mlabel = 'lcmv_contrast_threshold';
    case 2
        mlabel = 'lcmv_contrast_counting';
    case 3
        mlabel = 'lcmv_contrast_bootstrapping';
end

%%
save_dir = fullfile(data_save_dir,mlabel);

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    disp('Folder created successfully.');
else
    disp('Folder already exists.');
end

%% Compute LI (contrast)
close all
dtag = {'Ctrl';'Patn'};

for select_data = 1:length(dtag)
    
    disp(['Running ', mlabel, ' ', dtag{select_data}, '...'])
    
    cfg = [];
    cfg.sub_demog_data = sub_demog_data;
    cfg.select_data = select_data;
    S_data_sel = ecpfunc_select_data_contrast(cfg); 
    
    %- Subject-level LI (all 8 networks)
    cfg = [];
    cfg.S_data_sel = S_data_sel;
    cfg.BS_data_dir = BS_data_dir;
    cfg.Data_hcp_atlas = Data_hcp_atlas;
    cfg.idx_L = idx_L;
    cfg.idx_R = idx_R;
    cfg.thre = thre;
    cfg.wi = wi;
    cfg.data_save_dir = save_dir;
    cfg.method = mlabel;
    cfg.Threshtype = 3;
    [label_8net, LI_sub] = do_group_LI_net_lcmv_contrast(cfg);
end

%%
%- Plot LI subjects
cfg = [];
cfg.LI_sub = LI_sub;
cfg.wi = wi;
cfg.savefig = 1;
cfg.outdir = save_dir;
cfg.net_sel_mutiple_label = label_8net;
cfg.S_data_sel = S_data_sel;
cfg.network_sel = [1:3,6:10];
do_plot_group_lat(cfg);

cd(save_dir)
