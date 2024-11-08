clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/run')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/function')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/tools/helpful_tools/daviolinplot/daboxplot')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Main_scripts/run')
Run_setpath

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other/')

%%
flags = [];
flags.plot_atlasnetworks = 0;

%% HCP Atlas
src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim';
glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat';

glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';

cfg = struct('src_fname', src_fname, 'glass_dir', glass_dir, 'glass_atlas', glass_atlas, 'plotflag', 1);
Data_hcp_atlas = ecpfunc_hcp_atlas3(cfg);

%% SAVE HCP ATLAS NETWORKS
cfg = []; cfg.glass_dir = glass_dir; cfg.Data_hcp_atlas = Data_hcp_atlas; ecpfunc_hcp_atlas_save(cfg);

%% HCP Atlas, plotting
close all
cfg = [];
cfg.src_fname = src_fname;
cfg.network_sel = [1,2,6];
cfg.network_sel = [4];
cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.plotflag = 1;
cfg.fixedcolor = [0,0.7,0];
[idx_L, idx_R, src]  = do_plot_hcp_network(cfg);
net_rois = 'ftp';

disp(Data_hcp_atlas.groups_labels)


cfg = []; cfg.idx_L = idx_L; cfg.idx_R = idx_R; cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.src_fname = src_fname;

save_dir_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/atlas_roi';
outdir = save_dir_atlas;

cfg.export = 0; cfg.savedir = fullfile(outdir,'group');
cfg.network_sel = [1,2,6]; do_map_HCP_net_sel(cfg);
net_label = 'Fronto_tempro_pri';
cfg.network_sel = [11]; do_map_HCP_net_sel(cfg);



colorcode = {[0,114,189]; [217,83,25]; [237,177,32];[126,47,142]; [119,172,48]};
network_sel = [1,2,5,6,11]; % LI networks compared between MEG and fMRI

% Inspecting atlas networks common between MEG and fMRI
for i = 1:length(colorcode)
    cfg.network_sel = network_sel(i);
    cfg.fixedcolor = colorcode{i}/256;
    do_map_HCP_net_sel(cfg);
    title(Data_hcp_atlas.groups_labels{cfg.network_sel});
    
    cfg2 = [];
    cfg2.outdir = save_dir_atlas;
    filename = Data_hcp_atlas.groups_labels{cfg.network_sel};
    cfg2.filename = filename;
    cfg2.type = 'fig';
    %         do_export_fig(cfg2)
end