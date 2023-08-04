%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 06/27/2023

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
addpath('./run')
addpath('./function')
Run_setpath

%% Run Brainstorm
Run_BS

%%
Run_load_surface_template

%%
cfg = []; cfg.protocol = protocol;
S_data = ecpfunc_read_sourcemaps(cfg);

%% Subject demog details
cfg = []; cfg.subjs_3 = S_data.subjs_3; cfg.subjs_2 = S_data.subjs_2;
cfg.sFiles_3 = S_data.sFiles_3; cfg.sFiles_2 = S_data.sFiles_2;
sub_demog_data = ecpfunc_read_sub_demog(cfg);

%% TLE side (PT only)
patn_neuropsych_data = ecpfunc_read_patn_neuropsych();

%%
cfg = [];
cfg.sub_demog_data = sub_demog_data;
cfg.patn_neuropsych_data = patn_neuropsych_data;
sub_TLE_sub_data = ecpfunc_read_sub_TLE_sub(cfg);

%% Inter-subject (group) averaging,
disp('1: Anim, Ctrl')
disp('2: Anim, Patn')
disp('3: Symbol, Ctrl')
disp('4: Symbol, Patn')
select_data = input(':');

cfg = [];
cfg.sub_demog_data = sub_demog_data;
cfg.select_data = select_data;
S_data_sel = ecpfunc_select_data(cfg);

%% HCP Atlas
clc, close all
cfg = []; Data_hcp_atlas2 = ecpfunc_hcp_atlas2(cfg);
% cfg = []; Data_hcp_atlas = ecpfunc_hcp_atlas(cfg);

%% Time intervals (window)
cfg.strt = 0;
cfg.spt = 2;
% cfg.overlap = 0.01; cfg.linterval = 0.1;
cfg.overlap = 0.01;
cfg.linterval = 0.3;
wi  = do_time_intervals(cfg);

%%
thre = 0.5;

%%
src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/LI_subs';

%%
cfg = [];
cfg.src_fname = src_fname;
cfg.network_sel = [1,2,6];
% cfg.network_sel = [10];
cfg.Data_hcp_atlas = Data_hcp_atlas2;
[idx_L, idx_R, src]  = do_plot_hcp_network(cfg);
net_rois = 'ftp';

%% Whole-brain, avg, LI
close all
cfg = [];
cfg.idx_L = idx_L;
cfg.idx_R = idx_R;
cfg.thre = 0.8;
cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.S_data_sel = S_data_sel;
cfg.BS_data_dir = BS_data_dir;
cfg.wi = wi;
cfg.method = 'counting';
cfg.Threshtype = 3;
[idx_R_whole, idx_L_whole]  = do_LI_avg(cfg);
cfg.method = 'threshold';
cfg.thre = thre;
[idx_R_whole, idx_L_whole]  = do_LI_avg(cfg);
cfg.method = 'bootstrapping';
cfg.divs = 12; cfg.n_resampling = 2; cfg.RESAMPLE_RATIO = 0.8;
[idx_R_whole, idx_L_whole]  = do_LI_avg(cfg);

%% Network (all), avg, LI
% close all
cfg = [];
cfg.idx_L = idx_L;
cfg.idx_R = idx_R;
cfg.thre = thre;
cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.S_data_sel = S_data_sel;
cfg.BS_data_dir = BS_data_dir;
cfg.wi = wi;
cfg.thre = 0.8;
cfg.method = 'counting';
do_plot_LI_net_all(cfg)
cfg.thre = thre;
cfg.method = 'threshold';
do_plot_LI_net_all(cfg)
cfg.method = 'bootstrapping';
cfg.divs = 5; cfg.n_resampling = 2; 
cfg.RESAMPLE_RATIO = 1;
do_plot_LI_net_all(cfg)

%- Export the figure as a PDF/fig file
cfg = [];  cfg.outdir = fullfile(outdir,'group');  cfg.filename = [S_data_sel.s_tag, '-gmeanLI']; cfg.type = 'fig'; do_export_fig(cfg)
disp(fullfile(cfg.outdir,cfg.filename))

%% Network (all), avg, power - OPTIONAL
cfg = [];
cfg.idx_L = idx_L;
cfg.idx_R = idx_R;
cfg.thre = thre;
cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.S_data_sel = S_data_sel;
cfg.BS_data_dir = BS_data_dir;
cfg.wi = wi;
pow_lr = do_plot_pow_net_all(cfg);

%- Export figures as a PDF/fig file
cfg = [];
cfg.outdir = fullfile(outdir,'group');
cfg.filename = [S_data_sel.s_tag, '-pow'];
% cfg.type = 'pdf';
cfg.type = 'fig';
do_export_fig(cfg)

%%
disp(Data_hcp_atlas.groups_labels)

cfg = []; cfg.idx_L = idx_L; cfg.idx_R = idx_R; cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.src_fname = src_fname;

cfg.export = 1; cfg.savedir = fullfile(outdir,'group');
cfg.network_sel = [1,2,6]; do_map_HCP_net_sel(cfg);
net_label = 'Fronto_tempro_pri';

% Inspecting atlas areas.
for i=9:9
    cfg.network_sel = i; do_map_HCP_net_sel(cfg);title(Data_hcp_atlas.groups_labels{cfg.network_sel})
end

%% LI Subjects (network ROIs)
cfg = [];
cfg.idx_L_whole = idx_L_whole;
cfg.idx_R_whole = idx_R_whole;
cfg.thre = thre;
cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.S_data_sel = S_data_sel;
cfg.BS_data_dir = BS_data_dir;
cfg.wi = wi;
cfg.overwrite = 0;
cfg.data_save_dir = data_save_dir;
cfg.Threshtype = 3;
cfg.method = 'threshold';
[LI_sub, m_LI_sub, wi_sub_max] = do_sub_LI(cfg);
cfg.thre = 0.8;
cfg.method = 'counting';
[LI_sub, m_LI_sub, wi_sub_max] = do_sub_LI(cfg);
cfg.method = 'bootstrapping';
cfg.divs = 10;
[LI_sub, m_LI_sub, wi_sub_max] = do_sub_LI(cfg);

%% Export LI values
output_filename = fullfile(data_save_dir,'group', ['LI-', net_rois, '-', S_data_sel.s_tag, '.mat']);
save(output_filename,'m_LI_sub','wi_sub_max')

t1 = table(S_data_sel.sFiles_subid'); t1.Properties.VariableNames{'Var1'} = 'SubID';
t2 = table(m_LI_sub'); t2.Properties.VariableNames{'Var1'} = 'LI';
table_LI = [t1, t2];
writetable(table_LI, fullfile(data_save_dir,'group',['LI-', net_rois, '-', S_data_sel.s_tag,'.xlsx']));
disp(fullfile(data_save_dir,'group',['LI-', S_data_sel.s_tag,'.xlsx']))

sLI_sub = cell2mat(LI_sub');
figure, plot(mean(sLI_sub))
val = round(mean(wi(:,1),2),2);
set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
set(gca,'FontSize',8,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   400   1200   300]);
title('mean LI')

cfg = []; cfg.outdir = fullfile(outdir,'group');
cfg.filename = [S_data_sel.s_tag, '-mean_LI'];
cfg.type = 'fig';
do_export_fig(cfg)

cfg = [];
cfg.sub_sel = S_data_sel.sFiles_subid;
cfg.d_in = m_LI_sub;
cfg.tit = [S_data_sel.s_tag, '-sub-LI'];
do_barplot_LI(cfg)
set(gcf, 'Position', [1000   400   1000   300]);

cfg = []; cfg.outdir = fullfile(outdir,'group');
cfg.filename = [S_data_sel.s_tag, '-sub_LI'];
cfg.type = 'fig';
do_export_fig(cfg)

%% LI Subjects (Optimal ROIs & toi)
clc
cfg = [];
cfg.idx_L_whole = idx_L_whole;
cfg.idx_R_whole = idx_R_whole;
cfg.thre = thre;
cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.S_data_sel = S_data_sel;
cfg.BS_data_dir = BS_data_dir;
cfg.wi_sub_max = wi_sub_max;
do_sub_optimal_roi(cfg);

%% Inspecting: single-subject source subject maps (in time)
% cfg = [];
% cfg.S_data_sel = S_data_sel;
% cfg.wi_sub_max = wi_sub_max;
% cfg.BS_data_dir = BS_data_dir;
% cfg.src = src;
% do_sourcemap_time_sub(cfg);

%% time intervals of max LI
cfg = []; cfg.S_data_sel = S_data_sel; cfg.BS_data_dir = BS_data_dir;
cfg.wi_sub_max = wi_sub_max; cfg.src = src; do_sourcemap_time_sub_optimal_toi(cfg)

cfg = []; cfg.outdir = fullfile(outdir,'group');
cfg.filename = [S_data_sel.s_tag, '-max_LImap'];
% cfg.type = 'pdf';
cfg.type = 'fig';
do_export_fig(cfg)

% figure, bar(mean(wi_sub_max,2)), title('window')

%% mean intervals of max LI
tmp = mean(wi_sub_max);
cfg = []; cfg.S_data_sel = S_data_sel; cfg.BS_data_dir = BS_data_dir;
cfg.wi_sub_max = repmat(tmp, length(wi_sub_max), 1);
% cfg.wi_sub_max = repmat([250,650], length(wi_sub_max), 1);
cfg.src = src;
do_sourcemap_time_sub_optimal_toi(cfg)

cfg = []; cfg.outdir = fullfile(outdir,'group');
cfg.filename = [S_data_sel.s_tag, '-mean_LImap'];
% cfg.type = 'pdf';
cfg.type = 'fig';
do_export_fig(cfg)

%%
disp('1: threshold')
disp('2: counting')
disp('3: bootstrapping')
LI_method = input('LI_method sel:');
switch LI_method
    case 1
        mlabel = 'threshold';
    case 2
        mlabel = 'counting';
    case 3
        mlabel = 'bootstrapping';
end

%%
save_dir = fullfile(data_save_dir,'group_8net_300ms', mlabel);

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    disp('Folder created successfully.');
else
    disp('Folder already exists.');
end

%% Subject-level LI (all 8 networks)
cfg = [];
cfg.S_data_sel = S_data_sel;
cfg.BS_data_dir = BS_data_dir;
cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.idx_L = idx_L;
cfg.idx_R = idx_R;
cfg.thre = thre;
cfg.wi = wi;
cfg.data_save_dir = save_dir;
[label_8net, LI_sub, m_LI_max_sub, pow_sub] = do_group_LI_net(cfg);

%% Plot LI subjects
cfg = [];
cfg.net_sel_mutiple_label = label_8net;
cfg.LI_sub = LI_sub;
cfg.wi = wi;
cfg.subsel = 1;
cfg.plotflag = 1;
do_plot_sub_LI(cfg)

%% mean sub, all ROIs
cfg = [];
cfg.LI_sub = LI_sub;
cfg.wi = wi;
cfg.savefig = 1;
cfg.outdir = save_dir;
cfg.net_sel_mutiple_label = label_8net;
cfg.S_data_sel = S_data_sel;
cfg.network_sel = [1:3,6:8];
do_plot_group_lat(cfg);

cd(save_dir)


