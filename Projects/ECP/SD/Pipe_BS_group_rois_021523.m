
%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (preprocessing, source analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 01/11/2023

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
bs_path = '/opt/matlab_toolboxes/Brainstorm3_2021/brainstorm3'; % BS-2021
bs_path = '/opt/matlab_toolboxes/brainstorm3';

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
% if exist(fullfile(data_info_dir,'comments_subject.mat'),'file') == 2
%     load(fullfile(data_info_dir,'comments_subject.mat')),
% else
% clc
% close all,
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
save(fullfile(data_info_dir,'comments_subject.mat'),'Comment'),

%%
idx_anim = find(contains(Comment, 'Anim_')==1);
idx_symb = find(contains(Comment, 'Symbol_')==1);

sFiles_3 = sFiles_name(idx_anim);
sFiles_2 = sFiles_name(idx_symb);

%%
cd(data_info_dir)
% if exist(fullfile(data_info_dir,'subjects_ID.mat'),'file') == 2
%     load(fullfile(data_info_dir,'subjects_ID.mat')),
% else
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
% end

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
close all
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

idx_L = []; idx_R = []; idx_clr_L = []; idx_clr_R = [];
for i = 1:length(rois)
    idx = [];
    if find(strfind(rois{i}, 'L_')==1)
        idx_L = [idx_L; i];
    elseif find(strfind(rois{i}, 'R_')==1)
        idx_R = [idx_R; i];
    end
end
length(idx_R)
length(idx_L)

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
        sinput = sFiles_anim_hc;
    case 2
        sinput = sFiles_anim_pt;
    case 3
        sinput = sFiles_symb_hc;
    case 4
        sinput = sFiles_symb_pt;
end

%% Plot HCP atlas
cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.lat_index = [idx_L, idx_R];
cfg.rois = rois;
do_plot_HCP_atlas(cfg)

%% LI (network ROIs)
wi = []; w1 = 0.1; l = 0.1; ov = 0.01; j=1; %ov = l.*0.3
while w1+l < 2
    wi(j,:) = [w1, w1+l]; j=j+1; w1 = w1 + ov;
end
length(wi)
disp(wi)

%% Read Schaefer2018 atlas
% do_import_Schaefer

%% rois (from fMRI study)
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/glasser atlas';
% load(fullfile(glass_dir, 'LI_glasser_rois_48.mat'), 'glass_roi_name','glass_roi_R_idx','glass_roi_L_idx');
load(fullfile(glass_dir, 'LI_glasser_bilateral_rois_90.mat'), 'glass_roi_name','glass_roi_R_idx','glass_roi_L_idx');

cfg = [];
cfg.sinput = sFiles_anim_hc;
cfg.BS_data_dir = BS_data_dir; 
cfg.atlas = atlas;
cfg.thre = 0; 
cfg.wi = wi;
cfg.lat_index_L = glass_roi_L_idx; % glasser_lateral is not symmetric!
cfg.lat_index_R = glass_roi_R_idx; % glasser_lateral is not symmetric!
cfg.fplot = 1;
cfg.tit = 'test';
do_lat_analysis_asymetric(cfg);

%% network (from fMRI study)
% load(fullfile(glass_dir, 'LI_glasser_net_12.mat'));
load(fullfile(glass_dir, 'LI_glasser_bilateral_net_12.mat'));

disp('1: Angular'); disp('2: Frontal'); disp('3: Occipital');
disp('4: Other'); disp('5: PCingPrecun'); disp('6: Temporal');
network_ask  = input('enter network:');
switch network_ask
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
end

cfg = [];
cfg.sinput = sFiles_anim_hc;
cfg.BS_data_dir = BS_data_dir; 
cfg.atlas = atlas;
cfg.thre = 0; 
cfg.wi = wi;
cfg.lat_index_L = glass_net_L_idx{network_ask}; % glasser_lateral is not symmetric!
cfg.lat_index_R = glass_net_R_idx{network_ask}; % glasser_lateral is not symmetric!
cfg.fplot = 1;
cfg.tit = net_tag;
do_lat_analysis_asymetric(cfg);

groups_labels = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'};

clc
for i=1:length(groups_labels)
    disp([glass_net_L_label{i}, glass_net_R_label{i}])
end

cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
% cfg.lat_index_L = glass_net_L_name{network_ask};
% cfg.lat_index_R = glass_net_R_name{network_ask};
cfg.rois = rois;
cfg.group_labels = groups_labels;
cfg.group_members = glass_net_L_label;
cfg.roi_sel = network_ask;
do_plot_HCP6_atlas(cfg);
cfg.group_members = glass_net_R_label;
do_plot_HCP6_atlas(cfg);
% cfg.sel = 'whole'; % 'whole', 'left', 'right', 'roi';
% cfg.group_members = [glass_net_R_label,glass_net_L_label];
% do_plot_HCP6_atlas(cfg);

%%
[groups_labels, groups] = do_read_HCP_labels_bs;

cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.lat_index = [idx_L, idx_R];
cfg.rois = rois;
cfg.groups_labels = groups_labels;
cfg.groups = groups;
cfg.roi_sel = net_sel;
[idx23_L, idx23_R] = do_plot_HCP23_atlas(cfg);

%% Extract indecies
net_sel = [1,2,6];

cfg = [];
cfg.roi_network = 1:length(groups_labels); 
cfg.groups_labels = groups_labels; 
cfg.idx_L = idx_L; cfg.idx_R = idx_R;
cfg.net_sel = net_sel; 
[opt_idx_L,opt_idx_R] = do_network_indecies2(cfg);

%% LI Subjects (network ROIs)
cfg = [];
cfg.sinput = sFiles_anim_hc;
cfg.BS_data_dir = BS_data_dir;
cfg.atlas = atlas; cfg.thre = 0; cfg.fplot = 0; cfg.wi = wi;
cfg.lat_index = [opt_idx_L; opt_idx_R]';

ft_progress('init', 'text',     'please wait ...');
clear m_LI_sub LI_sub wi_sub_max
for i=1:length(sFiles_in)
    ft_progress(i/length(sFiles_in), 'Processing subjects %d from %d', i, length(sFiles_in));
    pause(0.1);
    cfg.sinput = sFiles_in{i};
    [LI,~, wi_max] = do_lat_analysis(cfg);
    LI_sub{i} = LI;
    m_LI_sub(i) = mean(LI);
    wi_sub_max(i,:) = wi_max;
end
ft_progress('close')

% Plot LIs
cfg = []; 
cfg.sub_sel = sub_sel; 
cfg.d_in = m_LI_sub;
cfg.tit = ['LIs (network): ', taskcon{J}, '-', subcon];
do_barplot_LI(cfg)

figure,
for i=1:length(LI_sub)
    plot(LI_sub{i}),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1500   500]);
    hold on
end

%%
disp('1: frontal')
disp('2: temporal')
disp('3: parietal')
disp('4: Fronto-temp-pari')
disp('5: Fronto-temp')
network_ask  = input('enter network:');
switch network_ask
    case 1
        net_sel = [12,21]; %[12,21:22];
        net_tag = 'frontal';
    case 2
        net_sel = [13,14,15];
        net_tag = 'temporal';
    case 3
        net_sel  = [16,17];
        net_tag = 'parietal';
    case 4
        % net_sel = [1:22];
        % net_sel = [12,13,14,15,16,17,21];
        % net_sel = [10,11,12,13,14,15,16,17,20, 21,22]
        % net_sel = [11,12,13,14,15,16,17,20, 21,22];
        net_sel = [11,12,13,14,16,17,20,21];
        net_tag = 'Fronto-temp-pari';
    case 5
        net_sel = [12,21,13,14,15];
        net_tag = 'Fronto-temp';
end

%%
% HCP 22 networks 
[groups_labels, groups] = do_read_HCP_labels_bs;

cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.lat_index = [idx_L, idx_R];
cfg.rois = rois;
cfg.groups_labels = groups_labels;
cfg.groups = groups;
cfg.roi_sel = net_sel;
[idx23_L, idx23_R] = do_plot_HCP23_atlas(cfg);

%%
net_sel_mutiple = {[12,21], [13,14,15], [16,17], [12,21,13,14,15], [11,12,13,14,16,17,20,21]};

cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.lat_index = [idx_L, idx_R];
cfg.rois = rois;
cfg.groups_labels = groups_labels;
cfg.groups = groups;
% cfg.roi_sel = [12,21]; do_plot_HCP23_atlas(cfg);
% cfg.roi_sel = [13,14,15]; do_plot_HCP23_atlas(cfg);
% cfg.roi_sel = [16,17]; do_plot_HCP23_atlas(cfg);
% cfg.roi_sel = [12,21,13,14,15]; do_plot_HCP23_atlas(cfg);
% cfg.roi_sel = [11,12,13,14,16,17,20,21]; do_plot_HCP23_atlas(cfg);
cfg.roi_sel = [1:22]; do_plot_HCP23_atlas(cfg);

% [idx_L32, idx_R32, groups_labels_num] = do_plot_HCP23_atlas(cfg);

%% Extract indecies
cfg = [];
cfg.roi_network = 1:length(groups_labels); 
cfg.groups_labels = groups_labels; 
cfg.idx_L = idx_L; cfg.idx_R = idx_R;
cfg.net_sel = net_sel; 
[opt_idx_L,opt_idx_R] = do_network_indecies(cfg);

%% Group avg. LI (network ROIs)
cfg = [];
cfg.sinput = sFiles_anim_hc;
cfg.BS_data_dir = BS_data_dir; 
cfg.atlas = atlas;
cfg.thre = 0; 
cfg.wi = wi;
cfg.lat_index = [opt_idx23_L; opt_idx_R]'; 
cfg.fplot = 1;
cfg.tit = net_tag;
[LI_hc_anim,rois_idx, L_max_hc_anim] = do_lat_analysis(cfg);

%% LI Subjects (network ROIs)
cfg = [];
cfg.sinput = sFiles_anim_hc;
cfg.BS_data_dir = BS_data_dir;
cfg.atlas = atlas; cfg.thre = 0; cfg.fplot = 0; cfg.wi = wi;
cfg.lat_index = [opt_idx23_L; opt_idx_R]';

ft_progress('init', 'text',     'please wait ...');
clear m_LI_sub LI_sub wi_sub_max
for i=1:length(sFiles_in)
    ft_progress(i/length(sFiles_in), 'Processing subjects %d from %d', i, length(sFiles_in));
    pause(0.1);
    cfg.sinput = sFiles_in{i};
    [LI,~, wi_max] = do_lat_analysis(cfg);
    LI_sub{i} = LI;
    m_LI_sub(i) = mean(LI);
    wi_sub_max(i,:) = wi_max;
end
ft_progress('close')

% Plot LIs
cfg = []; 
cfg.sub_sel = sub_sel; 
cfg.d_in = m_LI_sub;
cfg.tit = ['LIs (network): ', taskcon{J}, '-', subcon];
do_barplot_LI(cfg)

figure,
for i=1:length(LI_sub)
    plot(LI_sub{i}),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1500   500]);
    hold on
end

%% LI Subjects (Optimal ROIs & toi)
cfg = [];
cfg.sinput = sFiles_anim_hc;
cfg.BS_data_dir = BS_data_dir;
cfg.atlas = atlas; cfg.thre = 0; cfg.fplot = 0;
cfg.lat_index = [opt_idx_L; opt_idx_R]';

ft_progress('init', 'text',     'please wait ...');
clear m_LI_max_sub LI_sub 
for i=1:length(sFiles_in)
    ft_progress(i/length(sFiles_in), 'Processing subjects %d from %d', i, length(sFiles_in));
    pause(0.1);
    cfg.sinput = sFiles_in{i};
    cfg.wi = wi_sub_max(i,:);
    [LI,~, wi_max] = do_lat_analysis(cfg);
    LI_sub{i} = LI;
    m_LI_max_sub(i) = mean(LI);
%     wi_sub_max(i,:) = wi_max;
end
ft_progress('close')

% Plot LIs
cfg = []; 
cfg.sub_sel = sub_sel; 
cfg.d_in = m_LI_max_sub;
cfg.tit = ['LIs (network): ', taskcon{J}, '-', subcon];
do_barplot_LI(cfg)

% figure, bar(m_LI_sub, 0.4)
% set(gca,'Xtick', 1:length(sFiles_in),'XtickLabel',sub);
% set(gca,'FontSize',8,'XTickLabelRotation',90);
% set(gcf, 'Position', [1000   400   1000   500]);
% set(gca,'color','none');
% disp(wi_sub_max)

mean(wi_sub_max)
std(wi_sub_max)

%% Inspecting source subject maps (in time)
Run_sourcemap_time_sub

%%
Run_sourcemap_time_sub_optimal_toi

%% Subject-level LI
net_sel_mutiple_label = {'frontal', 'temporal', 'parietal', 'front-temporal(1)', 'front-temporal(2)'};
net_sel_mutiple = {[12,21], [13,14,15], [16,17], [12,21,13,14,15], [11,12,13,14,16,17,20,21]};

ft_progress('init', 'text',     'please wait ...');
clear m_LI_max_sub LI_sub
for j=1:length(net_sel_mutiple)
    ft_progress(j/length(net_sel_mutiple), 'Processing networks %d from %d', j, length(net_sel_mutiple));
    cfg = [];
    cfg.roi_network = roi_network; cfg.groups_labels = groups_labels;
    cfg.idx23_L = idx23_L; cfg.idx23_R = idx23_R;
    cfg.net_sel = net_sel_mutiple{j};
    [opt_idx_L,opt_idx_R] = do_network_indecies(cfg);
    
    %- Li Calc.
    cfg = [];
    cfg.sinput = sFiles_anim_hc;
    cfg.BS_data_dir = BS_data_dir;
    cfg.atlas = atlas; cfg.thre = 0; cfg.fplot = 0;
    cfg.lat_index = [opt_idx_L; opt_idx_R]';
    
    
    for i=1:length(sFiles_in)
        pause(0.1);
        cfg.sinput = sFiles_in{i};
        %     cfg.wi = wi_sub_max(i,:);
        cfg.wi = wi;
        [LI,~, wi_max] = do_lat_analysis(cfg);
        LI_sub(j,i,:) = LI;
        m_LI_max_sub(i) = mean(LI);
        %     wi_sub_max(i,:) = wi_max;
    end
    
end
ft_progress('close')

%%
% close all
for i=1:size(LI_sub,2)
    figure,
    for j=1:length(net_sel_mutiple)
        plot(squeeze(LI_sub(j,i,:))),
        val = round(mean(wi(:,1),2),2);
        set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
        set(gca,'FontSize',8,'XTickLabelRotation',90);
        set(gcf, 'Position', [1000   400   1100   300]);
        hold on
        legend(net_sel_mutiple_label)
        title(['sub:', num2str(i)])
    end
end

%%
mLI_sub = squeeze(mean(LI_sub,2));
figure,
for j=1:length(net_sel_mutiple)
    plot(mLI_sub(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   300]);
    hold on
    legend(net_sel_mutiple_label)
end
title('mean of LIs')

%% Group_mean-level LI
net_sel_mutiple_label = {'frontal', 'temporal', 'parietal', 'front-temporal', 'front-temporal-pari'};
net_sel_mutiple = {[12,21], [13,14,15], [16,17], [12,21,13,14,15], [11,12,13,14,16,17,20,21]};

clear m_LI_max_sub LI_sub
for j=1:length(net_sel_mutiple)
    
    cfg = [];
    cfg.roi_network = roi_network; cfg.groups_labels = groups_labels;
    cfg.idx23_L = idx23_L; cfg.idx23_R = idx23_R;
    cfg.net_sel = net_sel_mutiple{j};
    [opt_idx_L,opt_idx_R] = do_network_indecies(cfg);
    
    %- Li Calc.
    cfg = [];
    cfg.sinput = sFiles_anim_hc;
    cfg.BS_data_dir = BS_data_dir;
    cfg.atlas = atlas; cfg.thre = 0; cfg.fplot = 0;
    cfg.lat_index = [opt_idx_L; opt_idx_R]';
    
    cfg.sinput = sFiles_anim_hc;
    cfg.wi = wi;
    [LI,~, wi_max] = do_lat_analysis(cfg);
    LI_mean{j} = LI;    
end

figure,
for j=1:length(net_sel_mutiple)
    plot(LI_mean{j}),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1500   500]);
    hold on
    legend(net_sel_mutiple_label)
end

%% Response-time
clc,
rtdir = '/data/MEG/Research/aizadi/process';
d = rdir([rtdir,'/**/*RT.mat']);

RT_all = [];
for j=1:length(d)    
    RT_sel = load(d(j).name);
    RT_all.sub{j} = RT_sel.subj;
    RT_all.run{j} = RT_sel.run;
    RT_all.rt_time.animal(j) = mean(RT_sel.rt_time.animal);
    RT_all.rt_time.symbol(j) = mean(RT_sel.rt_time.symbol);
    RT_all.rt_time.both(j) = mean(RT_sel.rt_time.both);   
end

%%
[C,IA,IB] =  intersect(RT_all.sub,sub_sel);

%%
idx_val = [];
for i=1:length(C)
    idx = find(strcmp(C{i},sub_sel)==1);
    idx_val = [idx_val,idx];
end

idx_rt = [];
for i=1:length(C)
    idx = find(strcmp(C{i},RT_all.sub)==1);
    idx_rt = [idx_rt,idx];
end

%%
% rt_3 = RT_all.rt_time.animal(idx_rt);
% rt_2 = RT_all.rt_time.symbol(idx_rt);
% 
% crr_2_sel = crr_2(:,idx_val);
% crr_3_sel = crr_3(:,idx_val);
% 
% close all
% sub_crr_sel_2 = (crr_2_sel(:,1:2:end) + crr_2_sel(:,2:2:end))./2; figure, imagesc(sub_crr_sel_2'), title('2')
% sub_crr_sel_3 = (crr_3_sel(:,1:2:end) + crr_3_sel(:,2:2:end))./2; figure, imagesc(sub_crr_sel_3'), title('3')

%%
