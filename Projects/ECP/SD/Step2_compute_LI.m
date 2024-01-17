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

%%
flags = [];
flags.plot_atlasnetworks = 0;

%% HCP Atlas
clc, close all
src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/';
glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat';
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';

cfg = struct('src_fname', src_fname, 'glass_dir', glass_dir, 'glass_atlas', glass_atlas, 'plotflag', 1);
Data_hcp_atlas = ecpfunc_hcp_atlas2(cfg);

%% HCP Atlas, plotting
close all
cfg = [];
cfg.src_fname = src_fname;
cfg.network_sel = [1,2,6];
cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.plotflag = 1;
cfg.fixedcolor = [0,0.7,0];
[idx_L, idx_R, src]  = do_plot_hcp_network(cfg);
net_rois = 'ftp';

camlight(80,-10);
camlight(-80,-10);

%
disp(Data_hcp_atlas.groups_labels)


cfg = []; cfg.idx_L = idx_L; cfg.idx_R = idx_R; cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.src_fname = src_fname;

cfg.export = 1; cfg.savedir = fullfile(outdir,'group');
cfg.network_sel = [1,2,6]; do_map_HCP_net_sel(cfg);
net_label = 'Fronto_tempro_pri';

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

%% Run Brainstorm
Run_BS

%%
Run_load_surface_template

%%
LI_analysis_label = {'DICS_baseline','DICS_contrast','LCMV_basline','LCMV_contrast','DICS_anim'};

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
BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_all_subjects';

cfg = [];
cfg.protocol = protocol;
cfg.datadir = datadir;
cfg.BS_data_dir = BS_data_dir;

switch LI_analysis
    case {1,5}
%         cfg.datatag = 'wDICS_22_4_baseline';
        cfg.datatag = 'wDICS_baseline_18_4';
        S_data = ecpfunc_read_sourcemaps_dics(cfg);
    case 2
        cfg.datatag = 'wDICS_contrast_18_4';
%         cfg.datatag = 'wDICS_contrast_22_4';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 3
        S_data = ecpfunc_read_sourcemaps(cfg);
    case 4
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
    case {2, 4}
        cfg = []; cfg.subjs = S_data.subjs;
        cfg.sFiles = S_data.sFiles_32;
        sub_demog_data = ecpfunc_read_sub_demog_contrast(cfg);
end

%% Time intervals (window)
cfg.strt = 0;
cfg.spt = 2;
cfg.overlap = 0.01;
cfg.linterval = 0.3;
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
            [label_8net, LI_sub] = do_group_LI_net_baseline(cfg); % Anim vs. Baseline vs. Symb vs. Baseline
        end
        
    case {2, 4}
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
            [label_8net, LI_sub] = do_group_LI_net_contrast(cfg);
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
            
%             % Process: Average: Everything
%             bst_process('CallProcess', 'process_average', S_data_sel.sFiles_in', [], ...
%                 'avgtype',         1, ...  % Everything
%                 'avg_func',        1, ...  % Arithmetic average:  mean(x)
%                 'weighted',        0, ...
%                 'scalenormalized', 0);
            
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
            [label_8net, LI_sub] = do_group_LI_net_contrast(cfg); % Anim vs. Baseline vs. Symb vs. Baseline
        end
        
end

%%
clc, close all
% %- Plot LI subjects
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
