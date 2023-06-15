function Data_hcp_atlas = ecpfunc_hcp_atlas(~)

% ECP functions
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 05/31/2023

% atlas = load('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_corr.nii_362_updated.mat');
atlas = load('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat');
groups_labels = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'};

clear rois region region_lr
for i=1:length(atlas.Scouts)
    rois{i} = atlas.Scouts(i).Label;
    region{i} = atlas.Scouts(i).Region;
    region_lr{i} = atlas.Scouts(i).Region(1);
end

all_idx_L = []; all_idx_R = []; idx_clr_L = []; idx_clr_R = [];
for i = 1:length(rois)
    idx = [];
    if find(strfind(rois{i}, 'L_')==1)
        all_idx_L = [all_idx_L; i];
    elseif find(strfind(rois{i}, 'R_')==1)
        all_idx_R = [all_idx_R; i];
    end
end
% length(all_idx_R)
% length(all_idx_L)

%- Plot HCP atlas
% close all
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.index_L = all_idx_L;
cfg.index_R = all_idx_R;
cfg.rois = rois;
cfg.rois_sel = 1:180;
cfg.title = '';
do_plot_HCP_atlas(cfg)

%% Inspecting rois
%- Plot HCP atlas
% close all
% for i=[21]
%     disp(i)
%     cfg = [];
%     cfg.atlas = atlas;
%     cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
%     cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
%     cfg.index_L = all_idx_L;
%     cfg.index_R = all_idx_R;
%     cfg.rois = rois;
%     cfg.rois_sel = i;
%     cfg.title = num2str(i);
%     do_plot_HCP_atlas(cfg)
%     %     pause
% end

%% fixing frontal region (filling the missig regions)
k=1; clc
idx_sel_L = [];
idx_sel_roi_L = [];
for i=1:length(all_idx_L)
    tmp = atlas.Scouts(all_idx_L(i)).Region;
    if strcmp(tmp,'LF')
        idx_sel_L(k) = i;
        idx_sel_roi_L{k} = atlas.Scouts(all_idx_L(i)).Label;
        k=k+1;
    end
end

k=1; clc
idx_sel_R = [];
idx_sel_R_roi = [];
for i=1:length(all_idx_R)
    tmp = atlas.Scouts(all_idx_R(i)).Region;
    if strcmp(tmp,'RF')
        idx_sel_R(k) = i;
        idx_sel_R_roi{k} = atlas.Scouts(all_idx_R(i)).Label;
        k=k+1;
    end
end

% -inspection
% cfg = [];
% cfg.atlas = atlas;
% cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
% cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
% cfg.index_L = all_idx_L;
% cfg.index_R = all_idx_R;
% cfg.rois = rois;
% cfg.rois_sel = 166;
% cfg.title = '';
% do_plot_HCP_atlas(cfg)
% 
% idx_sel_roi_L{28}; % 164
% idx_sel_roi_L{30}; % 166
% 
% for i=idx_sel_L
%     disp(i)
%     cfg = [];
%     cfg.atlas = atlas;
%     cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
%     cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
%     cfg.index_L = all_idx_L;
%     cfg.index_R = all_idx_R;
%     cfg.rois = rois;
%     cfg.rois_sel = i;
%     cfg.title = num2str(i);
%     do_plot_HCP_atlas(cfg)
%     pause
% end

%% reading atlas rois (from fMRI study)
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/Glasser';
load(fullfile(glass_dir, 'LI_glasser_manual_net_12.mat'),'glass_net_L_label','glass_net_R_label');

%% updating frontal region
glass_net_L_label_update = glass_net_L_label;
glass_net_L_label_update{2} = [glass_net_L_label{2}; idx_sel_roi_L{28}; idx_sel_roi_L{30}];

glass_net_R_label_update = glass_net_R_label;
glass_net_R_label_update{2} = [glass_net_R_label{2}; idx_sel_R_roi{32}; idx_sel_R_roi{30}];


% net_sel = 2;
% 
% cfg = [];
% cfg.atlas = atlas;
% cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
% cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
% cfg.rois = rois;
% cfg.group_labels = groups_labels;
% cfg.group_members = glass_net_L_label_update;
% % Temporal
% cfg.roi_sel = net_sel;
% cfg.group_members = glass_net_L_label_update;
% cfg.group_members = glass_net_R_label_update;
% [idx_L, ~, ~] = do_plot_HCP6_atlas(cfg);

%% BTLA
net_sel = 6;
% 
% cfg = [];
% cfg.atlas = atlas;
% cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
% cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
% cfg.rois = rois;
% cfg.group_labels = groups_labels;
% cfg.group_members = glass_net_L_label_update;
% % Temporal
% cfg.roi_sel = net_sel;
% cfg.group_members = glass_net_L_label_update;
% [idx_L, ~, ~] = do_plot_HCP6_atlas(cfg);
% % cfg.group_members = glass_net_R_label;
% % [idx_R, ~, groups_labels_num] = do_plot_HCP6_atlas(cfg);
% 
% %- BTLA
% 
% rois_temp = idx_L{net_sel};

btla = [2,3,5,8,9,16,17,18,21,22];

% cfg = [];
% cfg.atlas = atlas;
% cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
% cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
% cfg.rois = rois;
% cfg.group_labels = groups_labels;
% glass_net_L_label_new = glass_net_L_label;
% glass_net_L_label_new{net_sel} = glass_net_L_label_new{net_sel}(btla);
% cfg.group_members = glass_net_L_label_new;
% cfg.roi_sel = net_sel;
% do_plot_HCP6_atlas(cfg);

BTLA_L_label = [];
BTLA_R_label = [];
for i=1:length(btla)
    BTLA_L_label = [BTLA_L_label; glass_net_L_label_update{net_sel}(btla(i))];
    BTLA_R_label = [BTLA_R_label; glass_net_R_label_update{net_sel}(btla(i))];
end

glass_net_L_label_update{7} = BTLA_L_label;
glass_net_R_label_update{7} = BTLA_R_label;

groups_labels{7} = 'BTLA';


%% Visual word form area VWFA
vw2 = [19,20]; net_sel = 6;
vw2 = [6,14,15,81]; net_sel = 4;

% cfg = [];
% cfg.atlas = atlas;
% cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
% cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
% cfg.rois = rois;
% cfg.group_labels = groups_labels;
% glass_net_L_label_new = glass_net_L_label;
% glass_net_L_label_new{net_sel} = glass_net_L_label_new{net_sel}(vw2);
% cfg.group_members = glass_net_L_label_new;
% cfg.roi_sel = net_sel;
% do_plot_HCP6_atlas(cfg);
% disp(glass_net_L_label_new)

VW_L_label = [];
VW_R_label = [];
for i=1:length(vw2)
    VW_L_label = [VW_L_label; glass_net_L_label_update{net_sel}(vw2(i))];
    VW_R_label = [VW_R_label; glass_net_R_label_update{net_sel}(vw2(i))];
end

glass_net_L_label_update{8} = VW_L_label;
glass_net_R_label_update{8} = VW_R_label;

groups_labels{8} = 'VWFA';

%%
k=1; clc
idx_sel_L = [];
idx_sel_roi_L = [];
for i=1:length(all_idx_L)
    tmp = atlas.Scouts(all_idx_L(i)).Region;
    if strcmp(tmp,'LT')
        idx_sel_L(k) = i;
        idx_sel_roi_L{k} = atlas.Scouts(all_idx_L(i)).Label;
        k=k+1;
    end
end

k=1; clc
idx_sel_R = [];
idx_sel_R_roi = [];
for i=1:length(all_idx_R)
    tmp = atlas.Scouts(all_idx_R(i)).Region;
    if strcmp(tmp,'RT')
        idx_sel_R(k) = i;
        idx_sel_R_roi{k} = atlas.Scouts(all_idx_R(i)).Label;
        k=k+1;
    end
end

% -inspection
% cfg = [];
% cfg.atlas = atlas;
% cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
% cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
% cfg.index_L = all_idx_L;
% cfg.index_R = all_idx_R;
% cfg.rois = rois;
% cfg.rois_sel = 1;
% cfg.title = '';
% do_plot_HCP_atlas(cfg)
% 
% idx_sel_roi_L{28}; % 164
% idx_sel_roi_L{30}; % 166
% 

%%
% for i=idx_sel_L
%     disp(i)
%     cfg = [];
%     cfg.atlas = atlas;
%     cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
%     cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
%     cfg.index_L = all_idx_L;
%     cfg.index_R = all_idx_R;
%     cfg.rois = rois;
%     cfg.rois_sel = i;
%     cfg.title = num2str(i);
%     do_plot_HCP_atlas(cfg)
%     pause
% end


%%
% for i=1:length(groups_labels)
%     net_sel = i;
%     cfg = [];
%     cfg.atlas = atlas;
%     cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
%     cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
%     cfg.rois = rois;
%     cfg.group_labels = groups_labels;
%     glass_net_L_label_new = glass_net_L_label_update;
%     glass_net_L_label_new{net_sel} = glass_net_L_label_new{net_sel};
%     cfg.group_members = glass_net_L_label_new;
%     cfg.roi_sel = net_sel;
%     do_plot_HCP6_atlas(cfg);
%     glass_net_R_label_new = glass_net_R_label_update;
%     glass_net_R_label_new{net_sel} = glass_net_R_label_new{net_sel};
%     cfg.group_members = glass_net_R_label_new;
%     do_plot_HCP6_atlas(cfg);
% end

%% BTLA
% Run_BTLA

%% Visual word form area VWFA
% Run_WVFA

%%
Data_hcp_atlas = [];
Data_hcp_atlas.glass_net_L_label = glass_net_L_label_update;
Data_hcp_atlas.glass_net_R_label = glass_net_R_label_update;
Data_hcp_atlas.groups_labels = groups_labels;
Data_hcp_atlas.atlas = atlas;
Data_hcp_atlas.rois = rois;

end