
%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (preprocessing, source analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 03/06/2023

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
cd('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD')
addpath('./run')
Run_setpath
addpath('./data')
addpath('./run')
% addpath(genpath('./functions'))
addpath(genpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions'))
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/Colormaps-from-MatPlotLib2.0')

results_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results';
addpath(ft_path); ft_defaults

%%
% adding BS path
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2021/brainstorm3'; % BS-2021
bs_path = '/opt/matlab_toolboxes/Brainstorm/brainstorm3_2019';

addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

BS_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/';
BS_data_dir = fullfile(BS_dir,'data_all_subjects');
protocol = fullfile(BS_dir, 'data_all_subjects/protocol.mat');

%%
% cd(results_dir)
src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
src = ft_read_headshape(src_fname);

src_S = [];
src_s.pnt = src.pos;
src_s.tri = src.tri;

%%
% db_reload_database('current',1)
load(protocol);
Subj_bs = ProtocolSubjects.Subject;

L = length(Subj_bs);
k = 1;
clear subjs_bs
for i=1:length(Subj_bs)
    if ~contains(Subj_bs(i).Name, 'Group_analysis')
        datafile{k} = Subj_bs(i).FileName;
        subjs_bs{k} = Subj_bs(i).Name;
        k=1+k;
    end
end
unq_bs_subj = unique(subjs_bs);

%%
no_anat = {'EC1036'
    'EC1037'
    'EC1038'
    'EC1040'
    'EC1045'
    'EC1049'
    'EC1061'
    'EC1065'
    'EC1085'
    'EC1092'
    'EC1094'
    'EC1096'
    'EC1110'
    'EC1111'
    'EC1112'
    'EC1141'
    'EC1153'
    'EC1162'
    'EC1090'};

sub_all1 = setdiff(unq_bs_subj,no_anat);

%%
data_info_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info';

cd(data_info_dir)
if exist(fullfile(data_info_dir,'comments_subject.mat'),'file') == 2
    load(fullfile(data_info_dir,'comments_subject.mat')),
else
    subj_del = [];
    
    cd(BS_data_dir)
    dd = rdir(fullfile('./Group_analysis/1_LCMV_Subjects/results_average*.mat'));
    for jj=1:length(dd), disp([num2str(jj),':',dd(jj).name]); end
    
    sFiles_name = [];
    for jj=1:length(dd)
        sFiles_name{jj} = fullfile(dd(jj).name(3:end));
    end
    Comment = []; k=1; kk=1;
    need_correction = [];
    for jj=1:length(sFiles_name)
        disp([num2str(jj), '/' , num2str(length(sFiles_name))])
        cd(BS_data_dir)
        tmp  = load(sFiles_name{jj});
        Comment{jj} = tmp.Comment;
        if length(tmp.Time) < 4000 && ~contains(Comment{jj}, 'Avg')
            disp(length(tmp.Time))
            need_correction(k) = jj;
            k=k+1;
        end
        if contains(Comment{jj}, 'Avg:')
            Comment_sel{kk} = sFiles_name{jj};
            kk=kk+1;
            disp(Comment{jj})
        end
    end
    
    removing_sFiles  = sFiles_name(need_correction);
    for i=1:length(removing_sFiles)
        %         removing_sFiles{i};
        [path, ~] = fileparts(removing_sFiles{i});
        if exist(path,'dir') == 7
            rmdir(fullfile(BS_data_dir, path),'s')
        end
    end
    
    idx_anim = find(contains(Comment, 'Anim_')==1);
    idx_symb = find(contains(Comment, 'Symbol_')==1);
    
    sFiles_3 = sFiles_name(idx_anim);
    sFiles_2 = sFiles_name(idx_symb);
    
    save(fullfile(data_info_dir,'comments_subject.mat'),'Comment','sFiles_name','sFiles_3','sFiles_2'),
end

%%
cd(data_info_dir)
if exist(fullfile(data_info_dir,'subjects_ID.mat'),'file') == 2
    load(fullfile(data_info_dir,'subjects_ID.mat')),
else
    cd(BS_data_dir)
    L = length(sFiles_3); k = 1;
    clear subjs_3
    for i=1:length(sFiles_3)
        %     disp([num2str(i), '/' , num2str(length(sFiles_3))])
        tmp = load(sFiles_3{i});
        idx = strfind(tmp.Comment,'_');
        subjs_3{k} = tmp.Comment(idx+1:idx+6);
        disp([subjs_3{k}, ': ', num2str(i), '/' , num2str(length(sFiles_3))])
        k=1+k;
    end
    unq_bs_subj_3 = unique(subjs_3);
    
    L = length(sFiles_2); k = 1;
    clear subjs_2
    for i=1:length(sFiles_2)
        %     disp([num2str(i), '/' , num2str(length(sFiles_3))])
        tmp = load(sFiles_2{i});
        idx = strfind(tmp.Comment,'_');
        subjs_2{k} = tmp.Comment(idx+1:idx+6);
        disp([subjs_2{k}, ': ', num2str(i), '/' , num2str(length(sFiles_2))])
        k=1+k;
    end
    unq_bs_subj_2 = unique(subjs_2);
    
    save(fullfile(data_info_dir,'subjects_ID.mat'),'subjs_2','subjs_3'),
end

%% Subject demog details
ECP_scriptdir = '/data/MEG/Research/ECP/Behavioural/processed';
load(fullfile(ECP_scriptdir,'sub_demog.mat'));

k=1; ib_3 = [];
for j=1:length(subjs_3)
    [~, ~,ib] = intersect(subjs_3{j},sub_demog_save(:,1));
    if ~isempty(ib)
        ib_3(k) = ib;
        k=k+1;
    end
end

k=1; ib_2 = [];
for j=1:length(subjs_2)
    [~, ~,ib] = intersect(subjs_2{j},sub_demog_save(:,1));
    if ~isempty(ib)
        ib_2(k) = ib;
        k=k+1;
    end
end

sub_cond_3 = sub_cond_val(ib_3); sub_cond_2 = sub_cond_val(ib_2);

idx_ctrl_3 = find(sub_cond_3 ==1); idx_patn_3 = find(sub_cond_3 ==2);
idx_ctrl_2 = find(sub_cond_2 ==1); idx_patn_2 = find(sub_cond_2 ==2);

sub_anim_hc = subjs_3(idx_ctrl_3);
sub_symb_hc = subjs_2(idx_ctrl_2);
sub_anim_pt = subjs_3(idx_patn_3);
sub_symb_pt = subjs_2(idx_patn_2);

sFiles_anim_patn = sFiles_3(idx_patn_3); sFiles_symb_patn = sFiles_2(idx_patn_2);
sFiles_anim_ctrl = sFiles_3(idx_ctrl_3); sFiles_symb_ctrl = sFiles_2(idx_ctrl_2);

%% Inter-subject (group) averaging,
clc
disp('1: Anim, Ctrl')
disp('2: Anim, Patn')
disp('3: Symbol, Ctrl')
disp('4: Symbol, Patn')
select_data = input(':');

switch select_data
    case 1
        J =1; subcon = 'Ctrl'; sFiles_in = sFiles_anim_ctrl; sub_sel = sub_anim_hc;  %- ctrl, Anim
    case 2
        J =1; subcon = 'Patn'; sFiles_in = sFiles_anim_patn; sub_sel = sub_anim_pt; %- Patn, Anim
    case 3
        J =2; subcon = 'Ctrl'; sFiles_in = sFiles_symb_ctrl; sub_sel = sub_symb_hc; %- Ctrl, Symb
    case 4
        J =2; subcon = 'Patn'; sFiles_in = sFiles_symb_patn; sub_sel = sub_symb_pt; %- Patn, Symb
end

taskcon = {'Anim', 'Symbol'};

%% HCP Atlas
% atlas = load('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_corr.nii_362_updated.mat');
atlas = load('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat');

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
length(all_idx_R)
length(all_idx_L)

%- Plot HCP atlas
close all
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.index_L = all_idx_L;
cfg.index_R = all_idx_R;
cfg.rois = rois;
cfg.rois_sel = 1:180;
cfg.title = '';
do_plot_HCP_atlas(cfg)

%%
%- Plot HCP atlas
close all
for i=[143,146,150,155,174]
    disp(i)
    cfg = [];
    cfg.atlas = atlas;
    cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
    cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
    cfg.index_L = all_idx_L;
    cfg.index_R = all_idx_R;
    cfg.rois = rois;
    cfg.rois_sel = i;
    cfg.title = num2str(i);
    do_plot_HCP_atlas(cfg)
    %     pause
end

%%
wi = []; w1 = 0; l = 0.1; ov = 0.01; j=1; %ov = l.*0.3
while w1+l < 2
    wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
end
length(wi)

%%
sFiles_anim_hc = 'Group_analysis/1_LCMV_Subjects/results_average_230110_2174.mat';
sFiles_anim_pt = 'Group_analysis/1_LCMV_Subjects/results_average_230110_2175.mat';
sFiles_symb_hc = 'Group_analysis/1_LCMV_Subjects/results_average_230111_1409.mat';
sFiles_symb_pt = 'Group_analysis/1_LCMV_Subjects/results_average_230113_1150.mat';
% sFiles_anim_hc_sFiles_anim_pt = 'Group_analysis/@intra/results_221226_1142.mat';

disp('1: Anim, Ctrl')
disp('2: Anim, Patn')
disp('3: Symbol, Ctrl')
disp('4: Symbol, Patn')
switch select_data
    case 1
        sinput = sFiles_anim_hc; s_tag = 'anim-hc';
    case 2
        sinput = sFiles_anim_pt; s_tag = 'anim-pt';
    case 3
        sinput = sFiles_symb_hc; s_tag = 'symb-hc';
    case 4
        sinput = sFiles_symb_pt; s_tag = 'symb-ppt';
end

%% rois (from fMRI study)
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/glasser atlas';
% load(fullfile(glass_dir, 'LI_glasser_bilateral_rois_90.mat'), 'glass_roi_name','glass_roi_R_idx','glass_roi_L_idx');
% load(fullfile(glass_dir, 'LI_glasser_bilateral_rois_90.mat'))

%% network (from fMRI study)
% load(fullfile(glass_dir, 'LI_glasser_net_12.mat'), 'glass_net_regions', 'glass_net_L_idx','glass_net_R_idx', 'glass_net_R_label', 'glass_net_L_label');
% load(fullfile(glass_dir, 'LI_glasser_bilateral_net_12.mat'));
load(fullfile(glass_dir, 'LI_glasser_manual_net_12.mat'));

%%
groups_labels = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'};

clc
for i=1:length(groups_labels)
    disp([glass_net_L_label{i}, glass_net_R_label{i}])
end

network_sel = [1,2,6];
% network_sel = [7,8];

cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.rois = rois;
cfg.group_labels = groups_labels;
cfg.group_members = glass_net_L_label;
cfg.roi_sel = network_sel;
do_plot_HCP6_atlas(cfg);
cfg.group_members = glass_net_R_label;
do_plot_HCP6_atlas(cfg);
% cfg.sel = 'whole'; % 'whole', 'left', 'right', 'roi';
% cfg.group_members = [glass_net_R_label,glass_net_L_label];
% do_plot_HCP6_atlas(cfg);

close all
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.rois = rois;
cfg.group_labels = groups_labels;
cfg.group_members = glass_net_L_label;
cfg.roi_sel = network_sel;
% net_sel = [6];
% cfg.roi_sel = net_sel;
cfg.group_members = glass_net_L_label;
[idx_L, ~, ~] = do_plot_HCP6_atlas(cfg);
cfg.group_members = glass_net_R_label;
[idx_R, ~, ~] = do_plot_HCP6_atlas(cfg);

%% Whole-brain, avg, LI
idx_R_whole = []; for i=1:length(idx_R), idx_R_whole = [idx_R_whole, idx_R{i}]; end
idx_L_whole = []; for i=1:length(idx_L), idx_L_whole = [idx_L_whole, idx_L{i}]; end

cfg = [];
cfg.sinput = sinput;
cfg.BS_data_dir = BS_data_dir;
cfg.atlas = atlas;
cfg.thre = 0;
cfg.wi = wi;
cfg.index_L = idx_L_whole; % glasser_lateral is not symmetric!
cfg.index_R = idx_R_whole; % glasser_lateral is not symmetric!
cfg.fplot = 1;
cfg.tit = 'Whole-brain, avg';%['Whole-brain, avg:', s_tag];
do_lat_analysis_asymetric(cfg);

%% BTLA
net_sel = [6];

close all
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.rois = rois;
cfg.group_labels = groups_labels;
cfg.group_members = glass_net_L_label;
% Temporal
cfg.roi_sel = net_sel;
cfg.group_members = glass_net_L_label;
[idx_L, ~, ~] = do_plot_HCP6_atlas(cfg);
cfg.group_members = glass_net_R_label;
[idx_R, ~, groups_labels_num] = do_plot_HCP6_atlas(cfg);

% %- Plot HCP atlas
rois_temp = idx_L{net_sel};
% close all
% for i=1:length(rois_temp)
%     %     disp(i)
%     %     cfg = [];
%     %     cfg.atlas = atlas;
%     %     cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
%     %     cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
%     %     cfg.index_L = all_idx_L;
%     %     cfg.index_R = all_idx_R;
%     %     cfg.rois = rois;
%     %     cfg.rois_sel = rois_temp(i);
%     %     cfg.title = num2str(rois_temp(i));
%     %     do_plot_HCP_atlas(cfg)
%
%     cfg = [];
%     cfg.atlas = atlas;
%     cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
%     cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
%     cfg.rois = rois;
%     cfg.group_labels = groups_labels;
%     glass_net_L_label_new = glass_net_L_label; glass_net_L_label_new{net_sel} = glass_net_L_label{net_sel}(i);
%     cfg.group_members = glass_net_L_label_new;
%     cfg.roi_sel = net_sel;
%     do_plot_HCP6_atlas(cfg);
% %     title(num2str(rois_temp(i)))
%     title(num2str((i)))
% %     glass_net_R_label_new = glass_net_R_label; glass_net_R_label_new{net_sel} = glass_net_R_label{net_sel}(i);
% %     cfg.group_members = glass_net_R_label_new;
% %     do_plot_HCP6_atlas(cfg);
% pause
% end

close all
btla = [2,3,5,8,9,16,17,18,21,22];
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.rois = rois;
cfg.group_labels = groups_labels;
glass_net_L_label_new = glass_net_L_label;
glass_net_L_label_new{net_sel} = glass_net_L_label_new{net_sel}(btla);
cfg.group_members = glass_net_L_label_new;
cfg.roi_sel = net_sel;
do_plot_HCP6_atlas(cfg);
title(num2str(rois_temp(i)))

BTLA_L_label = [];
BTLA_R_label = [];
for i=1:length(btla)
    BTLA_L_label = [BTLA_L_label; glass_net_L_label{net_sel}(btla(i))];
    BTLA_R_label = [BTLA_R_label; glass_net_R_label{net_sel}(btla(i))];
end

glass_net_L_label{7} = BTLA_L_label;
glass_net_R_label{7} = BTLA_R_label;

groups_labels{7} = 'BTLA';

close all
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.rois = rois;
cfg.group_labels = groups_labels;
cfg.group_members = glass_net_L_label;
% Temporal
cfg.roi_sel = 7;
cfg.group_members = glass_net_L_label;
[idx_L, ~, ~] = do_plot_HCP6_atlas(cfg);
cfg.group_members = glass_net_R_label;
[idx_R, ~, groups_labels_num] = do_plot_HCP6_atlas(cfg);


%% Visual word form area VWFA
net_sel = [4,6];

close all
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.rois = rois;
cfg.group_labels = groups_labels;
cfg.group_members = glass_net_L_label;
% Temporal
cfg.roi_sel = net_sel;
cfg.group_members = glass_net_L_label;
[idx_L, ~, ~] = do_plot_HCP6_atlas(cfg);
cfg.group_members = glass_net_R_label;
[idx_R, ~, groups_labels_num] = do_plot_HCP6_atlas(cfg);

%- Plot HCP atlas
% net_sel = 4;
% rois_other = idx_L{net_sel};
% close all
% for i=1:length(rois_other)
%     cfg = [];
%     cfg.atlas = atlas;
%     cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
%     cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
%     cfg.rois = rois;
%     cfg.group_labels = groups_labels;
%     glass_net_L_label_new = glass_net_L_label;
%     glass_net_L_label_new{net_sel} = glass_net_L_label{net_sel}(i);
%     cfg.group_members = glass_net_L_label_new;
%     cfg.roi_sel = net_sel;
%     do_plot_HCP6_atlas(cfg);
%     title(num2str((i)))
%     pause
% end

close all
net_sel = 4;
vw1 = [5,6,13,14,15,16,73,81,82,83,84,85];
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.rois = rois;
cfg.group_labels = groups_labels;
glass_net_L_label_new = glass_net_L_label;
glass_net_L_label_new{net_sel} = glass_net_L_label_new{net_sel}(vw1);
cfg.group_members = glass_net_L_label_new;
cfg.roi_sel = net_sel;
do_plot_HCP6_atlas(cfg);

%- Plot HCP atlas
% net_sel = 6;
% rois_temp = idx_L{net_sel};
% close all
% for i=1:length(rois_temp)
%     cfg = [];
%     cfg.atlas = atlas;
%     cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
%     cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
%     cfg.rois = rois;
%     cfg.group_labels = groups_labels;
%     glass_net_L_label_new = glass_net_L_label;
%     glass_net_L_label_new{net_sel} = glass_net_L_label{net_sel}(i);
%     cfg.group_members = glass_net_L_label_new;
%     cfg.roi_sel = net_sel;
%     do_plot_HCP6_atlas(cfg);
%     title(num2str((i)))
%     pause
% end

close all
vw2 = [19,20];
net_sel = 6;
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.rois = rois;
cfg.group_labels = groups_labels;
glass_net_L_label_new = glass_net_L_label;
glass_net_L_label_new{net_sel} = glass_net_L_label_new{net_sel}(vw2);
cfg.group_members = glass_net_L_label_new;
cfg.roi_sel = net_sel;
do_plot_HCP6_atlas(cfg);

VW_L_label = [];
VW_R_label = [];
for i=1:length(vw2)
    VW_L_label = [VW_L_label; glass_net_L_label{net_sel}(vw2(i))];
    VW_R_label = [VW_R_label; glass_net_R_label{net_sel}(vw2(i))];
end

glass_net_L_label{8} = VW_L_label;
glass_net_R_label{8} = VW_R_label;

groups_labels{8} = 'VW';

close all
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.rois = rois;
cfg.group_labels = groups_labels;
cfg.group_members = glass_net_L_label;
% Temporal
cfg.roi_sel = 8;
cfg.group_members = glass_net_L_label;
[idx_L, ~, ~] = do_plot_HCP6_atlas(cfg);
cfg.group_members = glass_net_R_label;
[idx_R, ~, groups_labels_num] = do_plot_HCP6_atlas(cfg);

%% Network, avg, LI
clc
disp('1: Angular'); disp('2: Frontal'); disp('3: Occipital');
disp('4: Other'); disp('5: PCingPrecun'); disp('6: Temporal');
disp('7: BTLA'); disp('8: Visual Word VWFA');
network_sel  = input('enter network:');
switch network_sel
    case 1
        net_tag = 'Angular';
    case 2
        net_tag = 'Frontal';
    case 3
        net_tag = 'Occipital';
    case 4
        net_tag = 'Other';
    case 5
        net_tag = 'PCingPrecun';
    case 6
        net_tag = 'Temporal';
    case 7
        net_tag = 'BTLA';
    case 8
        net_tag = 'VWFA';
end

cfg = [];
cfg.sinput = sinput;
cfg.BS_data_dir = BS_data_dir;
cfg.atlas = atlas;
cfg.thre = 0;
cfg.wi = wi;
cfg.index_L = idx_L{network_sel}; % glasser_lateral is not symmetric!
cfg.index_R = idx_R{network_sel}; % glasser_lateral is not symmetric!
cfg.fplot = 1;
cfg.tit = net_tag;
[LI,unq_roi_idx, LI_max] = do_lat_analysis_asymetric(cfg);

%%
close all
figure,
for i=1:length(idx_L)
    cfg = []; cfg.sinput = sinput; cfg.BS_data_dir = BS_data_dir;
    cfg.atlas = atlas; cfg.thre = 0; cfg.wi = wi;
    cfg.index_L = idx_L{i}; % glasser_lateral is not symmetric!
    cfg.index_R = idx_R{i}; % glasser_lateral is not symmetric!
    cfg.fplot = 0; cfg.tit = net_tag;
    [LI,unq_roi_idx, LI_max] = do_lat_analysis_asymetric(cfg);
    hold on
    plot(LI),
end
set(gcf, 'Position', [1000   400   1000   300]);
title(s_tag),
val = round(mean(wi(:,1),2),2);
legend(groups_labels)
set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
set(gca,'FontSize',8,'XTickLabelRotation',90);
xlabel('temporal windows (sec)')
ylabel('LI')
set(gca,'color','none');

%%
net_sel = [1,2,3,5,6,7,8];
net_sel = [8];

idx_R_all = [];
for i=1:length(net_sel)
    idx_R_all = [idx_R_all, idx_R{net_sel(i)}];
end

idx_L_all = [];
for i=1:length(net_sel)
    idx_L_all = [idx_L_all, idx_L{net_sel(i)}];
end

% close all
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'whole'; % 'whole', 'left', 'right', 'roi';
cfg.index_L = idx_L_all';
cfg.index_R = idx_R_all';
cfg.rois = rois;
cfg.rois_sel = 1:size(idx_L_all,2);
do_plot_HCP_atlas(cfg)

%% LI Subjects (network ROIs)
% cfg = [];
% cfg.sinput = sFiles_anim_hc;
% cfg.BS_data_dir = BS_data_dir;
% cfg.atlas = atlas;
% cfg.thre = 0;
% cfg.fplot = 0;
% cfg.wi = wi;
% cfg.index_L = idx_L_all;
% cfg.index_R = idx_R_all;
% 
% ft_progress('init', 'text',     'please wait ...');
% clear m_LI_sub LI_sub wi_sub_max
% for i=1:length(sFiles_in)
%     ft_progress(i/length(sFiles_in), 'Processing subjects %d from %d', i, length(sFiles_in));
%     pause(0.1);
%     cfg.sinput = sFiles_in{i};
%     [LI,~, wi_max] = do_lat_analysis_asymetric(cfg);
%     LI_sub{i} = LI;
%     m_LI_sub(i) = mean(LI);
%     wi_sub_max(i,:) = wi_max;
% end
% ft_progress('close')
% 
% % Plot LIs
% cfg = [];
% cfg.sub_sel = sub_sel;
% cfg.d_in = m_LI_sub;
% cfg.tit = ['LIs (network): ', taskcon{J}, '-', subcon];
% do_barplot_LI(cfg)
% set(gcf, 'Position', [1000   400   1000   300]);
% 
% 
% figure,
% for i=1:length(LI_sub)
%     plot(LI_sub{i}),
%     val = round(mean(wi(:,1),2),2);
%     set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
%     set(gca,'FontSize',8,'XTickLabelRotation',90);
%     %     set(gcf, 'Position', [1000   400   1500   500]);
%     hold on
% end
% set(gcf, 'Position', [1000   400   1000   300]);


%% LI Subjects (Optimal ROIs & toi)
% cfg = [];
% cfg.sinput = sFiles_anim_hc;
% cfg.BS_data_dir = BS_data_dir;
% cfg.atlas = atlas; cfg.thre = 0; cfg.fplot = 0;
% % cfg.lat_index = [opt_idx_L; opt_idx_R]';
% cfg.index_L = idx_L_all;
% cfg.index_R = idx_R_all;
% 
% ft_progress('init', 'text',     'please wait ...');
% clear m_LI_max_sub LI_sub
% for i=1:length(sFiles_in)
%     ft_progress(i/length(sFiles_in), 'Processing subjects %d from %d', i, length(sFiles_in));
%     pause(0.1);
%     cfg.sinput = sFiles_in{i};
%     cfg.wi = wi_sub_max(i,:);
%     [LI,~, wi_max] = do_lat_analysis_asymetric(cfg);
%     LI_sub{i} = LI;
%     m_LI_max_sub(i) = mean(LI);
%     %     wi_sub_max(i,:) = wi_max;
% end
% ft_progress('close')
% 
% % Plot LIs
% cfg = [];
% cfg.sub_sel = sub_sel;
% cfg.d_in = m_LI_max_sub;
% cfg.tit = ['LIs (network): ', taskcon{J}, '-', subcon];
% do_barplot_LI(cfg)
% set(gcf, 'Position', [1000   400   1000   300]);
% 
% % figure, bar(m_LI_sub, 0.4)
% % set(gca,'Xtick', 1:length(sFiles_in),'XtickLabel',sub);
% % set(gca,'FontSize',8,'XTickLabelRotation',90);
% % set(gcf, 'Position', [1000   400   1000   500]);
% % set(gca,'color','none');
% % disp(wi_sub_max)
% 
% mean(wi_sub_max)
% std(wi_sub_max)
% 
% %% Inspecting source subject maps (in time)
% Run_sourcemap_time_sub
% 
% %%
% Run_sourcemap_time_sub_optimal_toi

%% Subject-level LI
net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'BTLA'; 'VWFA'};
% net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'Whole'};
% idx_L{end+1} = idx_L_whole;idx_R{end+1} = idx_R_whole;

data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/LI_subs';

cd(data_save_dir)
savefilename = fullfile(data_save_dir,['LI_',s_tag, '.mat']);
if exist(savefilename,'file') == 2
    load(savefilename),
else
    
    ft_progress('init', 'text',     'please wait ...');
    clear m_LI_max_sub LI_sub
    for j=1:length(net_sel_mutiple_label)
        ft_progress(j/length(net_sel_mutiple_label), 'Processing networks %d from %d', j, length(net_sel_mutiple_label));
        
        %- Li Calc.
        cfg = [];
        cfg.sinput = sFiles_anim_hc;
        cfg.BS_data_dir = BS_data_dir;
        cfg.atlas = atlas; cfg.thre = 0; cfg.fplot = 0;
        cfg.index_L = idx_L{j};
        cfg.index_R = idx_R{j};
        
        
        for i=1:length(sFiles_in)
            pause(0.1);
            cfg.sinput = sFiles_in{i};
            %     cfg.wi = wi_sub_max(i,:);
            cfg.wi = wi;
            [LI,~, wi_max] = do_lat_analysis_asymetric(cfg);
            LI_sub(j,i,:) = LI;
            m_LI_max_sub(i) = mean(LI);
            %     wi_sub_max(i,:) = wi_max;
        end
        
    end
    ft_progress('close')
    save(savefilename,'Comment','LI_sub','m_LI_max_sub'),
end

%% LI subjects
% for i=1:size(LI_sub,2)
%     figure,
%     for j=1:length(net_sel_mutiple_label)
%         plot(squeeze(LI_sub(j,i,:))),
%         val = round(mean(wi(:,1),2),2);
%         set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
%         set(gca,'FontSize',8,'XTickLabelRotation',90);
%         set(gcf, 'Position', [1000   400   1100   300]);
%         hold on
%         title(['sub:', num2str(i)])
%     end
%     plot(squeeze(mean(LI_sub(:,i,:))),'LineWidth',3),
%     legend([net_sel_mutiple_label; 'mean'])
% end

%% mean sub, all ROIs
addpath('/data/MEG/Vahab/Github/MCW-MEGlab/FT/functions/External/brewermap');
network_sel = [1:3,6:8];
colr = distinguishable_colors(length(network_sel));

clc
close all
mLI_sub = squeeze(mean(LI_sub,2));
figure,
clear LI_val
for j=1:length(network_sel)
    LI_val(j,:) =  mLI_sub(network_sel(j),:);
    plot(LI_val(j,:), 'Color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   500]);
    hold on
end
p1 = plot(mean(LI_val),'LineWidth',3);
legend([net_sel_mutiple_label(network_sel); 'mean'])
title(['mean subject LIs: ', s_tag])
xlabel('time')
xlabel('time')

%%
clc
close all
mLI_sub = squeeze(mean(LI_sub,2));

colr_sub = distinguishable_colors(size(LI_sub,2));


clear std_dev
figure,
for j=1:length(network_sel)
    tmp = squeeze(LI_sub(network_sel(j),:,:));
    subplot(2,3,j)
    plot(tmp'),
    hold on
    mm(j,:) = mean(tmp);
    std_dev(j,:) = std(tmp);
    plot(mm(j,:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:5:length(wi),'XtickLabel',val(1:5:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gca,'color','none');
    title([s_tag, '-', net_sel_mutiple_label{network_sel(j)}])
    ylim([-0.6, 0.6])
%     xlim([0, 1.85])
    %     legend(net_sel_mutiple_label{j})
end

set(gcf, 'Position', [800   400   1600   500]);
% plot(mean(mLI_sub),'LineWidth',3),
% legend([net_sel_mutiple_label; 'mean'])

%% plot mean LI
figure,
for j=1:size(mm,1)
    hold on
    plot(mm(j,:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   500]);
end
lgnd = legend(net_sel_mutiple_label(network_sel));
title(['subject LIs: ', s_tag,])
set(gca,'color','none');
set(lgnd,'color','none');

%%
close all
figure,
for j=1:size(mm,1)
    a = mm(j,:);
    b = 1:numel(a);
    curve1 = a + std_dev(j,:);
    curve2 = a - std_dev(j,:);
    
    x2 = [b, fliplr(b)];
    inBetween = [curve1, fliplr(curve2)];
    hold on;
    h = plot(b, a, 'color', colr(j,:), 'LineWidth', 2);
    h2 = fill(x2, inBetween, colr(j,:), 'EdgeColor', 'none', 'facealpha', 0.1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off'; % make the legend for step plot off
    
end
lgnd = legend(net_sel_mutiple_label(network_sel));
set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
set(gca,'FontSize',8,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   400   1100   500]);
title(['subject LIs: ', s_tag,])
set(gca,'color','none');
set(lgnd,'color','none');
ylabel('LI')
xlabel('Time')


%%


