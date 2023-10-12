%% Logopenicppa dataset, Medical College of Wisconsin

% Script: BS Process (preprocessing, source analysis)
% Project: Logopenicppa_CRM
% Writtern by: Vahab YoussofZadeh
% Update: 05/28/2023

clear

addpath(genpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions'))

ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
addpath(ft_path);
ft_defaults

%% Reading data
BS_data_dir = '/data/MEG/Research/logopenicppa/BS_process/Logopenicppa_CRM_zscore';
cd(BS_data_dir)
dd = rdir(fullfile(BS_data_dir,'/data/Group_analysis/zscore/results*abs_ssmooth.mat'));

sFiles_name = []; for jj=1:length(dd), sFiles_name{jj} = fullfile(dd(jj).name); end


Comment = []; kk=1;
clear data_sub sess_num sub_ID sub_sess
for jj=1:length(sFiles_name)
    disp([num2str(jj), '/' , num2str(length(sFiles_name))])
    cd(BS_data_dir)
    tmp  = load(sFiles_name{jj});
    Comment{jj} = tmp.Comment;
    if contains(Comment{jj}, 'Avg:')
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
save(fullfile(data_info_dir,'comments_subject_zscore.mat'),'data_info'),

%% Atlas
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
disp('8: temp-pari')
disp('9: temp-pari (2)')
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
        net_sel = [17,15];
        net_tag = 'temporal 2';
    case 8
        net_sel = [13,14,15,16,17];
        %         net_sel = [13,14,15,17];
        net_tag = 'tempro-pari';
    case 9
        net_sel = [14, 15,16,17];
        net_tag = 'tempro-pari (2)';
end

[groups_labels, groups] = do_read_HCP_labels_bs;

cfg = [];
cfg.atlas = atlas;
% cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.lat_index = [idx_L, idx_R];
cfg.rois = rois;
cfg.groups_labels = groups_labels;
cfg.groups = groups;
cfg.roi_sel = net_sel;
[idx23_L, idx23_R] = do_plot_HCP23_atlas(cfg);

savedir = '/data/MEG/Research/logopenicppa/results/CRM/LI';
cfg = []; cfg.type = 'fig'; cfg.outdir = savedir; cfg.filename = ['Network_HCP_', net_tag]; do_export_fig(cfg)

%%
cfg = [];
cfg.roi_network = 1:length(groups_labels);
cfg.groups_labels = groups_labels;
cfg.idx_L = idx23_L; cfg.idx_R = idx23_R;
cfg.net_sel = net_sel;
[opt_idx_L,opt_idx_R] = do_network_indecies(cfg);

%% Time intervals (window)
clc
cfg = [];
cfg.strt = 0;
cfg.spt = 1;
cfg.overlap = 0.01;
cfg.linterval = 0.1;
wi  = do_time_intervals(cfg);

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
    cfg.fplot = 0;
    cfg.tit = net_tag;
    cfg.wi = wi;
    [LI,rois_idx, ~, pow] = do_lat_analysis(cfg);
    Li_sub(i,:) = LI;
    pow_sub{i} = pow;
end

%%
pow_sub_new = pow_sub;
pow_sub_new = pow_sub_new(2:end);

pow_sub_L = [];
pow_sub_R = [];
for i=1:length(pow_sub_new)
    pow_sub_L(i,:) = pow_sub_new{i}.left;
    pow_sub_R(i,:) = pow_sub_new{i}.right; 
end

LI_mean_left = [];
LI_mean_right = [];
tt = 20:60;
for i=1:length(pow_sub_new)
    mn_left = mean(pow_sub_L(i,tt));
    mn_right = mean(pow_sub_R(i,tt));
    pow_mean_left(i) = mn_left;
    pow_mean_right(i) = mn_right;
end

%%
% close all
sub_sess_new = sub_sess;
sub_sess_new(1) = [];

Li_sub_new = Li_sub;
Li_sub_new(1,:) = [];

data_info_new = data_info.sFiles_name;
data_info_new = data_info_new(2:end);

mLi_sub_new = mean(Li_sub_new(:,10:end),2);

% Plot LIs
cfg = [];
cfg.sub_sel = sub_sess_new;
cfg.d_in = mLi_sub_new;
cfg.tit = '';%['LIs (network): ', taskcon{J}, '-', subcon];
do_barplot_LI(cfg)
set(gcf, 'Position', [1000   500   1000   300]);

%%
wi_opt = [];
LI_max  =[];
LI_mean = [];
tonset = 1; % 300 ms
tt = 20:60;
t_sel = mean(wi,2);
t_sel1 = t_sel(tt);
disp(t_sel1)

for i=1:size(Li_sub_new,1)
    [mx, idx] = max(Li_sub_new(i,tt));
    mn = mean(Li_sub_new(i,tt));
    LI_max(i) = mx;
    LI_mean(i) = mn;
    tonset_wi = wi(tt,:);
    wi_opt(i,:) = tonset_wi(idx,:);
    [~, idx1] = min(abs(wi(:,1) - tonset_wi(idx)));
    wi_opt_idx(i) = idx1;
end
figure,plot(Li_sub_new'), size(Li_sub_new)

% figure,plot(Li_sub(1,:)),

%%
% close all
idx_subs = {1:4;5:8;9:12;13:16};

consistency = [];
for i=1:length(idx_subs)
    cr = corr(Li_sub_new(idx_subs{i},:)');
    consistency(i,:) = [mean(mean(cr(1:2))), mean(mean(cr(3:4)))];
end

consistency(1,:) = [consistency(1,2), consistency(1,1)];

figure, bar(consistency,0.2);
set(gca,'color','none'); axis square; xlabel('Subject'); ylabel('LI consistency'); box off

cfg = []; cfg.type = 'fig'; cfg.outdir = savedir; cfg.filename = ['LI_consist_subject_', net_tag]; do_export_fig(cfg)
cd(savedir)

disp(consistency)

%% LI analysis - sessions
savedir = '/data/MEG/Research/logopenicppa/results/CRM/LI';
colr = distinguishable_colors(4);
sub_sess_new2 = strrep(sub_sess_new, '_', '-');

for j=1:length(sub_sess_new2)
    sub_sess_new3{j} = sub_sess_new2{j}(end);
end

clc, close all
idx = 1:4; figure,subplot(4,3,1:2);
for i=1:length(idx), plot(mean(wi,2), Li_sub_new(idx(i),:)', 'color', colr(i,:)); hold on; end, box off
text(-0.2, 1.3, '(b)', 'Units', 'normalized', 'FontSize', 14);
% lgnd = legend(sub_sess_new3(idx),'Location', 'eastoutside'); set(lgnd,'color','none');
set(gca,'color','none');
tmp = LI_max(idx); tmp2 = [tmp(4)-tmp(3), tmp(2)-tmp(1)];
subplot(4,3,3), bar(100.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none');axis square, box off
text(-0.2, 1.3, '(c)', 'Units', 'normalized', 'FontSize', 14);
display(LI_max(idx));
ylim()

idx = 5:8; subplot(4,3,4:5);
for i=1:length(idx), plot(mean(wi,2), Li_sub_new(idx(i),:)', 'color', colr(i,:)); hold on; end, box off
% lgnd = legend(sub_sess_new3(idx),'Location', 'eastoutside'); set(lgnd,'color','none');
set(gca,'color','none');
tmp = LI_max(idx); tmp2 = [tmp(2)-tmp(1), tmp(4)-tmp(3)];
subplot(4,3,6), bar(100.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none'); axis square, box off
display(LI_max(idx))

idx = 9:12; subplot(4,3,7:8);
for i=1:length(idx), plot(mean(wi,2), Li_sub_new(idx(i),:)', 'color', colr(i,:)); hold on; end, box off
% lgnd = legend(sub_sess_new3(idx),'Location', 'eastoutside'); set(lgnd,'color','none');
set(gca,'color','none');
tmp = LI_max(idx); tmp2 = [tmp(2)-tmp(1), tmp(4)-tmp(3)];
subplot(4,3,9), bar(100.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none'); axis square, box off
display(LI_max(idx))

idx = 13:16; subplot(4,3,10:11);
for i=1:length(idx), plot(mean(wi,2), Li_sub_new(idx(i),:)', 'color', colr(i,:)); hold on; end, box off
xlabel('Time(sec)'); ylabel('Laterality')
lgnd = legend(sub_sess_new3(idx),'Location', 'northeast', 'Orientation', 'horizontal'); set(gca,'color','none'); set(lgnd,'color','none');
tmp = LI_max(idx); tmp2 = [tmp(2)-tmp(1), tmp(4)-tmp(3)];
subplot(4,3,12), bar(100.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]);
set(gca,'color','none'); axis square; xlabel('HD-tDCS'); ylabel({'mean laterality';'(% change from baseline)'}); box off
set(gcf, 'Position', [1000   400   600   500]);
display(LI_max(idx))
set(gca,'Xtick', 1:2,'XtickLabel',{'anodal ';'Sham'});

cfg = []; cfg.type = 'fig'; cfg.outdir = savedir; cfg.filename = ['LI_sessions_', net_tag]; do_export_fig(cfg)
cd(savedir)

%% LI
clc,
% close all
idx = 1:4;
figure,subplot(4,3,1:2);
tmp = Li_sub_new(idx(4),:) - Li_sub_new(idx(3),:); plot(t_sel1, tmp(tt), 'color', colr(1,:)); hold on;
tmp = Li_sub_new(idx(2),:) - Li_sub_new(idx(1),:); plot(t_sel1, tmp(tt), 'color', colr(2,:));
box off
text(-0.2, 1.3, '(b)', 'Units', 'normalized', 'FontSize', 14);
% lgnd = legend(sub_sess_new3(idx),'Location', 'eastoutside'); set(lgnd,'color','none');
set(gca,'color','none');
tmp = LI_mean(idx); tmp2 = [tmp(4)-tmp(3), tmp(2)-tmp(1)];
LI_anodal_sham(1) = -diff(tmp2);
subplot(4,3,3), bar(100.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none');axis square, box off
text(-0.2, 1.3, '(c)', 'Units', 'normalized', 'FontSize', 14);
display(100.*tmp2);

idx = 5:8; subplot(4,3,4:5);
tmp = Li_sub_new(idx(2),:) - Li_sub_new(idx(1),:); plot(t_sel1, tmp(tt), 'color', colr(1,:)); hold on;
tmp = Li_sub_new(idx(4),:) - Li_sub_new(idx(3),:); plot(t_sel1, tmp(tt), 'color', colr(2,:));
box off
set(gca,'color','none');
tmp = LI_mean(idx); tmp2 = [tmp(2)-tmp(1), tmp(4)-tmp(3)];
LI_anodal_sham(2) = -diff(tmp2);
subplot(4,3,6), bar(100.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none'); axis square, box off
display(100.*tmp2);

idx = 9:12; subplot(4,3,7:8);
tmp = Li_sub_new(idx(2),:) - Li_sub_new(idx(1),:); plot(t_sel1, tmp(tt), 'color', colr(1,:)); hold on;
tmp = Li_sub_new(idx(4),:) - Li_sub_new(idx(3),:); plot(t_sel1, tmp(tt), 'color', colr(2,:));
set(gca,'color','none');
box off
tmp = LI_mean(idx); tmp2 = [tmp(2)-tmp(1), tmp(4)-tmp(3)];
LI_anodal_sham(3) = -diff(tmp2);
subplot(4,3,9), bar(100.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none'); axis square, box off
display(100.*tmp2);

idx = 13:16; subplot(4,3,10:11);
tmp = Li_sub_new(idx(2),:) - Li_sub_new(idx(1),:); plot(t_sel1, tmp(tt), 'color', colr(1,:)); hold on;
tmp = Li_sub_new(idx(4),:) - Li_sub_new(idx(3),:); plot(t_sel1, tmp(tt), 'color', colr(2,:));
box off
xlabel('Time(sec)'); ylabel('Laterality')
lgnd = legend(sub_sess_new3(idx),'Location', 'northeast', 'Orientation', 'horizontal'); set(gca,'color','none'); set(lgnd,'color','none');
tmp = LI_mean(idx); tmp2 = [tmp(2)-tmp(1), tmp(4)-tmp(3)];
LI_anodal_sham(4) = -diff(tmp2);
subplot(4,3,12), bar(100.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]);
set(gca,'color','none'); axis square; xlabel('HD-tDCS'); ylabel({'mean laterality';'(% change from baseline)'}); box off
set(gcf, 'Position', [1000   400   600   500]);
display(100.*tmp2);
set(gca,'Xtick', 1:2,'XtickLabel',{'anodal ';'Sham'});

cfg = []; cfg.type = 'fig'; cfg.outdir = savedir; cfg.filename = ['LI_sessions_', net_tag]; do_export_fig(cfg)
cd(savedir)
disp(cfg.filename)

%%
figure, bar(100.*LI_anodal_sham,0.2);
set(gca,'color','none'); axis square; xlabel('Subject'); ylabel('LI anodal - sham'); box off

disp(100.*LI_anodal_sham)

cfg = []; cfg.type = 'fig'; cfg.outdir = savedir; cfg.filename = ['LI_diff_', net_tag]; do_export_fig(cfg)
cd(savedir)

%% Power - Left
clc,
close all
idx = 1:4;
figure,subplot(4,3,1:2);
tmp = pow_sub_L(idx(4),:) - pow_sub_L(idx(3),:); plot(t_sel1, tmp(tt), 'color', colr(1,:)); hold on;
tmp = pow_sub_L(idx(2),:) - pow_sub_L(idx(1),:); plot(t_sel1, tmp(tt), 'color', colr(2,:));
box off
text(-0.2, 1.3, '(b)', 'Units', 'normalized', 'FontSize', 14);
% lgnd = legend(sub_sess_new3(idx),'Location', 'eastoutside'); set(lgnd,'color','none');
set(gca,'color','none');
tmp = pow_mean_left(idx); tmp2 = [tmp(4)-tmp(3), tmp(2)-tmp(1)];
pow_anodal_sham(1) = -diff(tmp2);
subplot(4,3,3), bar(1000.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none');axis square, box off
text(-0.2, 1.3, '(c)', 'Units', 'normalized', 'FontSize', 14);
display(1000.*tmp2);

idx = 5:8; subplot(4,3,4:5);
tmp = pow_sub_L(idx(2),:) - pow_sub_L(idx(1),:); plot(t_sel1, tmp(tt), 'color', colr(1,:)); hold on;
tmp = pow_sub_L(idx(4),:) - pow_sub_L(idx(3),:); plot(t_sel1, tmp(tt), 'color', colr(2,:));
box off
set(gca,'color','none');
tmp = pow_mean_left(idx); tmp2 = [tmp(2)-tmp(1), tmp(4)-tmp(3)];
pow_anodal_sham(2) = -diff(tmp2);
subplot(4,3,6), bar(1000.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none'); axis square, box off
display(1000.*tmp2);

idx = 9:12; subplot(4,3,7:8);
tmp = pow_sub_L(idx(2),:) - pow_sub_L(idx(1),:); plot(t_sel1, tmp(tt), 'color', colr(1,:)); hold on;
tmp = pow_sub_L(idx(4),:) - pow_sub_L(idx(3),:); plot(t_sel1, tmp(tt), 'color', colr(2,:));
set(gca,'color','none');
box off
tmp = pow_mean_left(idx); tmp2 = [tmp(2)-tmp(1), tmp(4)-tmp(3)];
pow_anodal_sham(3) = -diff(tmp2);
subplot(4,3,9), bar(1000.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]); set(gca,'color','none'); axis square, box off
display(1000.*tmp2);

idx = 13:16; subplot(4,3,10:11);
tmp = pow_sub_L(idx(2),:) - pow_sub_L(idx(1),:); plot(t_sel1, tmp(tt), 'color', colr(1,:)); hold on;
tmp = pow_sub_L(idx(4),:) - pow_sub_L(idx(3),:); plot(t_sel1, tmp(tt), 'color', colr(2,:));
box off
xlabel('Time(sec)'); ylabel('Power')
lgnd = legend(sub_sess_new3(idx),'Location', 'northeast', 'Orientation', 'horizontal'); set(gca,'color','none'); set(lgnd,'color','none');
tmp = pow_mean_left(idx); tmp2 = [tmp(2)-tmp(1), tmp(4)-tmp(3)];
pow_anodal_sham(4) = -diff(tmp2);
subplot(4,3,12), bar(1000.*tmp2, 0.2, 'FaceColor',[.5 .5 .5]);
set(gca,'color','none'); axis square; xlabel('HD-tDCS'); ylabel({'mean power';'(% change from baseline)'}); box off
set(gcf, 'Position', [1000   400   600   500]);
display(1000.*tmp2);
set(gca,'Xtick', 1:2,'XtickLabel',{'anodal ';'Sham'});

cfg = []; cfg.type = 'fig'; cfg.outdir = savedir; cfg.filename = ['Pow_L_sessions_', net_tag]; do_export_fig(cfg)
cd(savedir)
disp(cfg.filename)

%%
figure, bar(1000.*pow_anodal_sham,0.2);
set(gca,'color','none'); axis square; xlabel('Subject'); ylabel('Pow anodal - sham'); box off

disp(1000.*pow_anodal_sham)

cfg = []; cfg.type = 'fig'; cfg.outdir = savedir; cfg.filename = ['Pow_L_diff_', net_tag]; do_export_fig(cfg)
cd(savedir)

%%
cfg = [];
cfg.sub_sel = sub_sess_new;
cfg.d_in = LI_max;
cfg.tit = '';%['LIs (network): ', taskcon{J}, '-', subcon];
do_barplot_LI(cfg)
set(gcf, 'Position', [1000   500   1000   300]);

%% Export maps with optimal LI
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3'; %BS-2021
addpath(bs_path);
brainstorm
disp('choose DB from BS, then enter!');
pause

%%
disp(wi_opt)
BSpath = '/data/MEG/Research/logopenicppa/BS_process/Logopenicppa_CRM_zscore/data/';
for i=1:length(data_info_new)
    cfg_in = [];
    cfg_in.side_sel = 7;
    cfg_in.BSpath = BSpath;
    cfg_in.svdir = '/data/MEG/Research/logopenicppa/results/CRM/LI_zscore';
    idx = strfind(data_info_new{i}, 'Group_analysis');
    cfg_in.fname = data_info_new{i}(idx:end);
    cfg_in.sname = [net_tag, '-', sub_sess_new{i}];
    cfg_in.seltime = mean(wi_opt(i,:));
    cfg_in.imgres = 1; % high-res = 1, low-res = 2;
    disp(wi_opt(i))
    do_export_images_BS(cfg_in),
end

%%