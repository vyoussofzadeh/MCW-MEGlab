%% Logopenicppa dataset, Medical College of Wisconsin

% Script: BS Process (preprocessing, source analysis)
% Project: Logopenicppa_CRM
% Writtern by: Vahab YoussofZadeh
% Update: 04/25/2023

addpath(genpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions'))

ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
% ft_path = '/opt/matlab_toolboxes/ft_packages/fieldtrip_20190419'
addpath(ft_path);
ft_defaults

%% Reading data
BS_data_dir = '/data/MEG/Research/logopenicppa/BS_process/Logopenicppa_CRM_new';
cd(BS_data_dir)
% '/data/MEG/Research/logopenicppa/BS_process/Logopenicppa_CRM_new/data/32037_002_1/@intra';
% dd1 = rdir(fullfile(BS_data_dir,'/data/Group_analysis/@intra/results_average*ssmooth.mat'));
% dd = rdir(fullfile(BS_data_dir,'/data/Group_analysis/@intra/results_average*.mat'));
dd = rdir(fullfile(BS_data_dir,'/data/Group_analysis/smoothed_avg_300_1000ms/results*.mat'));

% list1 = {dd.name};
% list2 = {dd1.name};
% 
% dd2 = setdiff(list1, list2);

% for jj=1:length(dd2), disp([num2str(jj),':',dd2{jj}]); end
% 
% sFiles_name = [];
% for jj=1:length(dd2), sFiles_name{jj} = fullfile(dd2{jj}); end

for jj=1:length(dd), disp([num2str(jj),':',dd(jj).name]); end

sFiles_name = [];
for jj=1:length(dd), sFiles_name{jj} = fullfile(dd(jj).name); end

Comment = []; kk=1;
clear data_sub sess_num sub_ID sub_sess
for jj=1:length(sFiles_name)
    disp([num2str(jj), '/' , num2str(length(sFiles_name))])
    cd(BS_data_dir)
    tmp  = load(sFiles_name{jj});
    Comment{jj} = tmp.Comment;
    if contains(Comment{jj}, '/Avg:')
        Comment_sel{kk} = sFiles_name{jj};
        disp(Comment{jj});
        tkz = tokenize(Comment{jj},'/');
        data_sub{kk} = tkz{1};
        tkz1 = tokenize(tkz{1},'_');
        sub_ID{kk} = tkz1{2};
        sess_num{kk} = tkz1{3};
        sub_sess{kk} = [sub_ID{kk},'_',sess_num{kk}];
        kk=kk+1;
    end
end

data_info_dir = '/data/MEG/Research/logopenicppa/Scripts';

data_info = [];
data_info.Comment = Comment;
data_info.sub_sess = sub_sess;
data_info.data_sub = data_sub;
data_info.sub_ID = sub_ID;
data_info.sess_num = sess_num;
data_info.sFiles_name = sFiles_name;

disp('----')
disp(data_info)
save(fullfile(data_info_dir,'comments_subject_avg.mat'),'data_info'),

%% read atlas
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
close all
disp('1: frontal')
disp('2: temporal')
disp('3: parietal')
disp('4: Fronto-temp-pari')
disp('5: Fronto-temp')
disp('6: all')
disp('7: temporal 2')
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
        net_sel = [12,13,14,15, 21];
        net_tag = 'Fronto-temp';
    case 6
        net_sel = 1:22;
        net_tag = 'Whole-brain';
    case 7
        net_sel = [14,15];
        net_tag = 'temporal 2';
end

%
[groups_labels, groups] = do_read_HCP_labels_bs;

cfg = [];
cfg.atlas = atlas;
% cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.lat_index = [idx_L, idx_R];
cfg.rois = rois;
cfg.groups_labels = groups_labels;
cfg.groups = groups;
cfg.roi_sel = net_sel;
[idx23_L, idx23_R] = do_plot_HCP23_atlas(cfg);

savedir = '/data/MEG/Research/logopenicppa/results/CRM/LI';
cfg = []; cfg.type = 'fig'; cfg.outdir = savedir; cfg.filename = 'Network_HCP'; do_export_fig(cfg)

%%
cfg = [];
cfg.roi_network = 1:length(groups_labels);
cfg.groups_labels = groups_labels;
cfg.idx_L = idx23_L; cfg.idx_R = idx23_R;
cfg.net_sel = net_sel;
[opt_idx_L,opt_idx_R] = do_network_indecies(cfg);

%% Time intervals (window)
% cfg = [];
% cfg.strt = 0;
% cfg.spt = 1;
% cfg.overlap = 0.01;
% wi  = do_time_intervals(cfg);

%% LI Subjects (network ROIs)
clear Li_sub
ft_progress('init', 'text',     'please wait ...');
for i=1:length(data_info.sFiles_name)
    ft_progress(i/length(data_info.sFiles_name), 'Processing subjects %d from %d', i, length(data_info.sFiles_name));
    
    cfg = [];
    cfg.sinput = data_info.sFiles_name{i};
    cfg.BS_data_dir = '';
    cfg.atlas = atlas;
    cfg.thre = 0.5;
    cfg.lat_index = [opt_idx_L; opt_idx_R]'; % network_based
    %     cfg.lat_index = [idx_L, idx_R]; % Whole-brain (voxel_based)
    cfg.fplot = 0;
    cfg.tit = net_tag;
    cfg.wi = [1 1];
    %     [LI,rois_idx]
    cfg.fplot = 0;
    [LI,unq_roi_idx, LI_max, pow] = do_lat_analysis(cfg);
    Li_sub(i,:) = LI;
end

%%
% close all
sub_sess_new = sub_sess;
% sub_sess_new(1) = [];

Li_sub_new = Li_sub;
% Li_sub_new(1,:) = [];

data_info_new = data_info.sFiles_name;
% data_info_new = data_info_new(2:end);

% mLi_sub_new = mean(Li_sub_new(:,10:end),2);

% Plot LIs
cfg = [];
cfg.sub_sel = sub_sess_new;
cfg.d_in = Li_sub_new;
cfg.tit = '';%['LIs (network): ', taskcon{J}, '-', subcon];
do_barplot_LI(cfg)
set(gcf, 'Position', [1000   500   1000   300]);

%%
% wi_opt = [];
% tonset = 1; % 300 ms
% for i=1:size(Li_sub_new,1)
%     [mx, idx] = max(Li_sub_new(i,tonset:end));
%     LI_max(i) = mx;
%     tonset_wi = wi(tonset:end,:);
%     wi_opt(i) = tonset_wi(idx);
%     [~, idx1] = min(abs(wi(:,1) - tonset_wi(idx)));
%     wi_opt_idx(i) = idx1;
% end
% figure,plot(Li_sub_new'), size(Li_sub_new)
% 
% % figure,plot(Li_sub(1,:)),

%% LI analysis - sessions

% savedir = '/data/MEG/Research/logopenicppa/results/CRM/LI';
% colr = distinguishable_colors(4);
% sub_sess_new2 = strrep(sub_sess_new, '_', '-');
% 
% clc, close all
% idx = 1:4; figure,subplot(4,3,1:2);
% for i=1:length(idx), plot(mean(wi,2), Li_sub_new(i,:)', 'color', colr(i,:)); hold on; end, box off
% text(-0.2, 1.3, '(a)', 'Units', 'normalized', 'FontSize', 14);
% lgnd = legend(sub_sess_new2(idx),'Location', 'eastoutside'); set(gca,'color','none'); set(lgnd,'color','none');
% subplot(4,3,3), bar(LI_max(idx), 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none');axis square, box off
% text(-0.2, 1.3, '(b)', 'Units', 'normalized', 'FontSize', 14);
% 
% idx = 5:8; subplot(4,3,4:5); 
% for i=1:length(idx), plot(mean(wi,2), Li_sub_new(i,:)', 'color', colr(i,:)); hold on; end, box off
% lgnd = legend(sub_sess_new2(idx),'Location', 'eastoutside'); set(gca,'color','none'); set(lgnd,'color','none');
% subplot(4,3,6), bar(LI_max(idx), 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none'); axis square, box off
% 
% idx = 9:12; subplot(4,3,7:8); 
% for i=1:length(idx), plot(mean(wi,2), Li_sub_new(i,:)', 'color', colr(i,:)); hold on; end, box off
% lgnd = legend(sub_sess_new2(idx),'Location', 'eastoutside'); set(gca,'color','none'); set(lgnd,'color','none');
% subplot(4,3,9), bar(LI_max(idx), 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none'); axis square, box off
% 
% idx = 13:16; subplot(4,3,10:11); 
% for i=1:length(idx), plot(mean(wi,2), Li_sub_new(i,:)', 'color', colr(i,:)); hold on; end, box off
% xlabel('Time(sec)'); ylabel('Laterality') 
% lgnd = legend(sub_sess_new2(idx),'Location', 'eastoutside'); set(gca,'color','none'); set(lgnd,'color','none');
% subplot(4,3,12), bar(LI_max(idx), 0.2, 'FaceColor',[.5 .5 .5]); 
% set(gca,'color','none'); axis square; xlabel('Session'); ylabel('mean laterality'); box off
% set(gcf, 'Position', [1000   400   600   500]);
% 
% cfg = []; cfg.type = 'fig'; cfg.outdir = savedir; cfg.filename = 'LI_sessions'; do_export_fig(cfg)

%%
% cfg = [];
% cfg.sub_sel = sub_sess_new;
% cfg.d_in = LI_max;
% cfg.tit = '';%['LIs (network): ', taskcon{J}, '-', subcon];
% do_barplot_LI(cfg)
% set(gcf, 'Position', [1000   500   1000   300]);

%% Export maps with optimal LI 
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3'; %BS-2021
addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

%%
disp(wi_opt')
BSpath = '/data/MEG/Research/logopenicppa/BS_process/Logopenicppa_CRM_new/data/';
for i=1:length(data_info_new)
    cfg_in = [];
    cfg_in.side_sel = 2;
    cfg_in.BSpath = BSpath;
    cfg_in.svdir = '/data/MEG/Research/logopenicppa/results/CRM/STG';
    idx = strfind(data_info_new{i}, 'Group_analysis');
    cfg_in.fname = data_info_new{i}(idx:end);
    disp(wi_opt(i))
    do_export_images_BS(cfg_in),
end

%%
% 0.1600
% 0.2800
% 0.3400

% 0.5400

% 0.4500

% 0.1400

% 0.1500

% 0.3000

% 0.4900

% 0.8200

% 0.3000

% 0.3900

% 0.3400

% 0.1900

% 0.2600

% 0.4200

%%
% 1: 0.8000
% 2: 0.3500
% 3: 0.4500
% 4: 0.5300
% 5: 0.5300
% 6: 0.0900
% 7: 0.4300
% 8: 0.5800
% 9: 0.1500
% 10: 0.3700
% 11: 0.7200
% 12: 0.6900
% 13: 0.4300
% 14: 0.2200
% 15: 0.3000
% 16: 0.3300

%%
% cfg_in = [];
% cfg_in.side_sel = 2;
% cfg_in.svdir = '/data/MEG/Research/logopenicppa/results';
% cfg_in.BSpath = '/data/MEG/Research/logopenicppa/BS_process/Logopenicppa_CRM/data';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1728.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1723_32037_001_1.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1712_32037_002_1.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1713_32037_002_2.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1714_32037_002_3.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1714_32037_002_4.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1715_32037_004_1.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1716_32037_004_2.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1716_32037_004_3.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1717_32037_004_4.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1718_32037_005_1.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1718_32037_005_2.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1719_32037_005_3.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1720_32037_005_4.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1720_32037_006_1.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1721_32037_006_2.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_220317_1722_32037_006_3.mat';
% cfg_in.fname = 'Group_analysis/@intra/results_average_230220_1619_32037_006_4.mat';
% do_export_images_BS(cfg_in),

%%
% figure,plot(Li_sub(1,:)),
% val = round(mean(wi(:,1),2),2);
% set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
% set(gca,'FontSize',8,'XTickLabelRotation',90);
% set(gcf, 'Position', [1000   400   1500   500]);
% title([cfg.tit, ' - ', tmp.Comment]),
% xlabel('temporal windows (sec)')
% ylabel('LI')
% set(gca,'color','none');
% 
% %%
% close all
% cfg = [];
% cfg.LI_sub = Li_sub;
% cfg.wi = wi;
% cfg.savefig = 1;
% cfg.outdir = fullfile(outdir,'group'); 
% cfg.net_sel_mutiple_label = net_sel_mutiple_label;
% cfg.S_data_sel = S_data_sel;
% do_plot_group_lat(cfg);

