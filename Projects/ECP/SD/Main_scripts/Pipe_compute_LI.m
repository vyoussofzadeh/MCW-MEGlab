%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 06/27/2023

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/run')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/function')
Run_setpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other/')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/tools/helpful_tools/daviolinplot/daboxplot');

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

%% SAVE ATLAS NETWORKS
cfg = []; cfg.glass_dir = glass_dir; cfg.Data_hcp_atlas = Data_hcp_atlas; ecpfunc_hcp_atlas_save(cfg);

%%
% Lateral network
% close all
cfg = [];
cfg.src_fname = src_fname;
cfg.network_sel = [11];
cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.plotflag = 0;
cfg.fixedcolor = [0,0.7,0];
[idx_L, idx_R, src]  = do_plot_hcp_network(cfg);
net_rois = 'ftp';
disp(Data_hcp_atlas.groups_labels)

%
cfg = []; cfg.idx_L = idx_L; cfg.idx_R = idx_R; cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.src_fname = src_fname;

% cfg.export = 1; cfg.savedir = fullfile(outdir,'group');
% cfg.network_sel = [1,2,6]; do_map_HCP_net_sel(cfg);
% net_label = 'Fronto_tempro_pri';

save_dir_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/atlas_roi';

colorcode = {[0,114,189]; [217,83,25]; [237,177,32];[126,47,142]; [119,172,48]};
network_sel = [1,2,5,6,11]; % LI networks compared between MEG and fMRI

if flags.plot_atlasnetworks == 1
    
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
        do_export_fig(cfg2)
    end
    
    
    % Inspecting all atlas networks areas (MEG only)
    for i = 1:11
        cfg.network_sel = i;
        cfg.fixedcolor = colorcode/256;
        do_map_HCP_net_sel(cfg);
        title(Data_hcp_atlas.groups_labels{cfg.network_sel});
        
        cfg2 = [];
        cfg2.outdir = save_dir_atlas;
        filename = Data_hcp_atlas.groups_labels{cfg.network_sel};
        cfg2.filename = filename;
        cfg2.type = 'fig';
        do_export_fig(cfg2)
    end
end

%% BS
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3';

addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

BS_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/';
protocol = fullfile(BS_dir, 'data_full/protocol.mat');

%%
Run_load_surface_template

%%
% LI_analysis_label = {'DICS_baseline','DICS_contrast','LCMV_baseline','LCMV_contrast','DICS_anim', 'DICS_contrast_prestim', 'dSPM_contrast'};
LI_analysis_label = {'DICS_indirect','DICS_directcontrast','LCMV_anim_vs_Symb','-','DICS_anim', 'DICS_contrast_prestim', 'dSPM_contrast'};

for i = 1:length(LI_analysis_label)
    disp([num2str(i) ') ' LI_analysis_label{i}]);
end

LI_analysis = input('');

%%
LI_method_label = {'Magnitude', 'Counting','Bootstrapping'};

disp('1: Magnitude')
disp('2: Counting')
disp('3: Bootstrapping')
LI_method = input(':');

%%
datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';
% BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_full';
BS_data_dir = fullfile(BS_dir,'data_full');

cfg = [];
cfg.protocol = protocol;
cfg.datadir = datadir;
cfg.BS_data_dir = BS_data_dir;

switch LI_analysis
    case 1
        cfg.datatag = 'wDICS_baseline_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics(cfg);
    case 2
        cfg.datatag = 'wDICS_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 3
        cfg.datamask = fullfile('./Group_analysis/LCMV/results_average*.mat');
        S_data = ecpfunc_read_sourcemaps(cfg);
        %     case 4
        %         cfg.datamask = fullfile('./Group_analysis/LCMV/results_abs*.mat');
        %         S_data = ecpfunc_read_sourcemaps_contrast(cfg);
    case 5
        protocol = fullfile(BS_dir, 'data/protocol.mat');
        BS_data_dir = fullfile(BS_dir,'data');
        cfg.protocol = protocol;
        cfg.BS_data_dir = BS_data_dir;
        cfg.datatag = 'wDICS_baseline_18_4';
        S_data = ecpfunc_read_sourcemaps_dics_anim(cfg);
    case 6
        cfg.datatag = 'wDICS_contrast_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 7
        cfg.datatag = 'dSPM_contrast';
        S_data = ecpfunc_read_sourcemaps_contrast(cfg);
end

%% TLE side (PT only)
patn_neuropsych_data = ecpfunc_read_patn_neuropsych();

%% Subject demog details
switch LI_analysis
    case {1,3 ,5}
        cfg = []; cfg.subjs_3 = S_data.subjs_3; cfg.subjs_2 = S_data.subjs_2;
        cfg.sFiles_3 = S_data.sFiles_3; cfg.sFiles_2 = S_data.sFiles_2;
        sub_demog_data = ecpfunc_read_sub_demog(cfg);
    case {2, 6, 7}
        cfg = []; cfg.subjs = S_data.subjs;
        cfg.sFiles = S_data.sFiles_32;
        sub_demog_data = ecpfunc_read_sub_demog_contrast(cfg);
end

%% Time intervals (window)
cfg =[];
cfg.strt = -0.5;
cfg.spt = 2;
cfg.overlap = 0.05;
% cfg.overlap = 0.01;
% cfg.overlap = 0.15;
cfg.linterval = 0.3;
% cfg.linterval = 0.1;
wi  = do_time_intervals(cfg);

%%
thre = 0.5;

%%
save_dir = fullfile(data_save_dir,LI_analysis_label{LI_analysis}, LI_method_label{LI_method});
checkOrCreateDir(save_dir)
cd(save_dir)

%% Compute LI
mlabel = [LI_analysis_label{LI_analysis}, '_', LI_method_label{LI_method}];

switch LI_analysis
    case {1,3}
        dtag = {'Anim, Ctrl';'Anim, Patn'; 'Symbol, Ctrl'; 'Symbol, Patn'};
        dtag_val = {[1,3]; [2,4]};
        
        for select_data = 1:length(dtag_val)
            
            disp(['Running ', mlabel, ' ', dtag{dtag_val{select_data}}, '...'])
            
            cfg = [];
            cfg.sub_demog_data = sub_demog_data;
            cfg.select_data = dtag_val{select_data}(1);
            S_data_sel_A = ecpfunc_select_data(cfg);
            cfg.select_data = dtag_val{select_data}(2);
            S_data_sel_B = ecpfunc_select_data(cfg);
            
            [C,IA,IB] = intersect(S_data_sel_A.sFiles_subid, S_data_sel_B.sFiles_subid);
            
            S_data_sel_A1 = S_data_sel_A;
            S_data_sel_A1.sFiles_subid = S_data_sel_A1.sFiles_subid(IA);
            S_data_sel_A1.sFiles_in = S_data_sel_A1.sFiles_in(IA);
            
            S_data_sel_B1 = S_data_sel_B;
            S_data_sel_B1.sFiles_subid = S_data_sel_B1.sFiles_subid(IB);
            S_data_sel_B1.sFiles_in = S_data_sel_B1.sFiles_in(IB);
            
            S_data_sel = {S_data_sel_A1; S_data_sel_B1};
            
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
            cfg.doavg = 1;
            [label_8net, LI_sub] = do_group_LI_net_baseline(cfg); % Anim vs. Baseline vs. Symb vs. Baseline
        end
        
    case {2, 6, 7}
        dtag = {'Ctrl';'Patn'};
        
        for select_data = 1:length(dtag)
            
            disp(['Running ', mlabel, ' ', dtag{select_data}, '...'])
            
            cfg = [];
            cfg.sub_demog_data = sub_demog_data;
            cfg.select_data = select_data;
            S_data_sel = ecpfunc_select_data_contrast(cfg);
            
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
            cfg.doavg = 0;
            [label_8net, LI_sub] = do_group_LI_net_contrast_dics(cfg);
        end
    case 5
        dtag = {'Anim, Ctrl';'Anim, Patn'};
        dtag_val = 1:2;
   
        for select_data = 1:length(dtag_val)
            
            disp(['Running ', mlabel, ' ', dtag{dtag_val(select_data)}, '...'])
            
            cfg = [];
            cfg.sub_demog_data = sub_demog_data;
            cfg.select_data = dtag_val(select_data);
            S_data_sel = ecpfunc_select_data(cfg);
            
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
            cfg.doavg = 0;
            cfg.parcellaion = 0;
            [label_8net, LI_sub] = do_group_LI_net_contrast(cfg); % Anim vs. Baseline vs. Symb vs. Baseline
        end    
end

%% Plot LI subjects
% clc, close all
cfg = [];
cfg.LI_sub = LI_sub;
cfg.wi = wi;
cfg.savefig = 0;
cfg.outdir = save_dir;
cfg.net_sel_mutiple_label = label_8net;
cfg.S_data_sel = S_data_sel;
cfg.network_sel = [11];
do_plot_group_lat(cfg);
%
% cd(save_dir)
