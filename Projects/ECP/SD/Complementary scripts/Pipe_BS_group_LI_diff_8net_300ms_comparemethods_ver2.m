%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 05/31/2023

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

cfg = [];
cfg.sub_demog_data = sub_demog_data;
cfg.select_data = 1; S_data_anim_hc = ecpfunc_select_data(cfg);
cfg.select_data = 2; S_data_anim_pt = ecpfunc_select_data(cfg);
cfg.select_data = 3; S_data_symb_hc = ecpfunc_select_data(cfg);
cfg.select_data = 4; S_data_symb_pt = ecpfunc_select_data(cfg);

%%
cfg = []; cfg.strt = 0; cfg.spt = 2; cfg.overlap = 0.01; cfg.linterval = 0.3;
wi  = do_time_intervals(cfg);

%%
net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'BTLA'; 'VWFA'};

%%
network_sel = [1:3,6:8];
colr = distinguishable_colors(length(network_sel));

%%
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/LI_subs/group_8net_300ms';
cd(data_save_dir)

%%
close all
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
%
cd(fullfile(data_save_dir, mlabel))

LI_anim_hc = load('LI_anim-hc');
LI_anim_pt = load('LI_anim-pt');

LI_symb_hc = load('LI_symb-hc');
LI_symb_pt = load('LI_symb-pt');

%%
clc
d_in =[];

for LI_method=1:3
    switch LI_method
        case 1
            mlabel = 'threshold';
        case 2
            mlabel = 'counting';
        case 3
            mlabel = 'bootstrapping';
    end
    d_label = 'LI_anim-hc'; d_in.(mlabel).LI_anim_hc = load(fullfile(data_save_dir, mlabel,[d_label, '.mat']));
    d_label = 'LI_anim-pt'; d_in.(mlabel).LI_anim_pt = load(fullfile(data_save_dir, mlabel,[d_label, '.mat']));
    d_label = 'LI_symb-hc'; d_in.(mlabel).LI_symb_hc = load(fullfile(data_save_dir, mlabel,[d_label, '.mat']));
    d_label = 'LI_symb-pt'; d_in.(mlabel).LI_symb_pt = load(fullfile(data_save_dir, mlabel,[d_label, '.mat']));   
    
end

% cfg = [];
% cfg.d_in = d_in;
% cfg.network_sel = network_sel;
% do_compare_methods(cfg)

%%
close all

pow = [];
if LI_method == 1
    for j=1:size(LI_anim_hc.pow_sub,1)
        for i=1:size(LI_anim_hc.pow_sub,2)
            pow.left_anim_hc(j,i,:) = LI_anim_hc.pow_sub(j,i).left;
            pow.right_anim_hc(j,i,:) = LI_anim_hc.pow_sub(j,i).right;
            pow.left_symb_hc(j,i,:) = LI_symb_hc.pow_sub(j,i).left;
            pow.right_symb_hc(j,i,:) = LI_symb_hc.pow_sub(j,i).right;
        end
    end
end

mPow_sub1 = squeeze(mean(pow.left_anim_hc,2));
mPow_sub2 = squeeze(mean(pow.left_symb_hc,2)); tag = [mlabel,'; anim vs. symb, hc, left'];

%-
clc
mPow_sub_hc = mPow_sub1 - mPow_sub2;

figure,
clear LI_val
for j=1:length(network_sel)
    hold on
    do_createPlot(mPow_sub_hc(network_sel(j),:), val, colr(j,:), net_sel_mutiple_label(network_sel), [mlabel,'; hc - pt'], 'LI')
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(tag)
ylabel('Pow')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');


mPow_sub1 = squeeze(mean(pow.right_anim_hc,2));
mPow_sub2 = squeeze(mean(pow.right_symb_hc,2)); tag = [mlabel,'; anim vs. symb, hc, right'];

%-
clc
mPow_sub_hc = mPow_sub1 - mPow_sub2;

figure,
clear LI_val
for j=1:length(network_sel)
    hold on
    do_createPlot(mPow_sub_hc(network_sel(j),:), val, colr(j,:), net_sel_mutiple_label(network_sel), [mlabel,'; hc - pt'], 'LI')
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(tag)
ylabel('Pow')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');

%% Anim HC
% close all
mLI_sub1 = squeeze(mean(LI_anim_hc.LI_sub,2)); tag = [mlabel, '; anim hc'];

clc, mLI_sub_hc = mLI_sub1;

figure,
clear LI_val
for j=1:length(network_sel)
    hold on
    plot(mLI_sub_hc(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   500]);
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(tag)
ylabel('LI')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');

%% anim vs. symb - HC
% close all
mLI_sub1 = squeeze(mean(LI_anim_hc.LI_sub,2));
mLI_sub2 = squeeze(mean(LI_symb_hc.LI_sub,2)); tag = [mlabel,'; anim vs. symb, hc'];

clc
mLI_sub_hc = mLI_sub1 - mLI_sub2;

figure,
% subplot(3,1,1)
clear LI_val
for j=1:length(network_sel)
    hold on
    plot(mLI_sub_hc(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   500]);
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(tag)
ylabel('LI')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');

%%
[sub_pt,IA,IB] = intersect(S_data_anim_pt.sFiles_subid, S_data_symb_pt.sFiles_subid);

LI_anim_pt_val = LI_anim_pt.LI_sub(:,IA,:);
LI_symb_pt_val = LI_symb_pt.LI_sub(:,IB,:);


mLI_sub1 = squeeze(mean(LI_anim_pt_val,2));
mLI_sub2 = squeeze(mean(LI_symb_pt_val,2)); 

tag = 'anim vs. symb, pt';

mLI_sub_pt = mLI_sub1 - mLI_sub2;

figure,
% subplot(3,1,2)
clear LI_val
for j=1:length(network_sel)
    hold on
    plot(mLI_sub_pt(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   500]);
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(tag)
ylabel('LI')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');

%%
mLI_sub_diff = mLI_sub_hc - mLI_sub_pt; tag = [mlabel,'; hc - pt'];

figure,
% subplot(3,1,3)
clear LI_val
for j=1:length(network_sel)
    hold on
    plot(mLI_sub_diff(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   500]);
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(tag)
ylabel('LI')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');

set(gcf, 'Position', [1000   400   1000   900]);

%%
% mLI_sub_mdiff = mean(mLI_sub_hc,1) - mean(mLI_sub_pt,1); tag = 'hc - pt, mean';
% 
% figure,
% plot(mLI_sub_mdiff,'LineWidth',3, 'color','k'),
% val = round(mean(wi(:,1),2),2);
% set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
% set(gca,'FontSize',8,'XTickLabelRotation',90);
% set(gcf, 'Position', [1000   400   1100   300]);
% title(tag)
% set(gca,'color','none');
% set(lgnd,'color','none');

%%
patn_neuropsych_tle = ecpfunc_read_patn_neuropsych_tle();
TLESide = patn_neuropsych_tle.TLESide; SUBNO = patn_neuropsych_tle.SUBNO;

for i=1:length(sub_pt) 
    SUBNO_anim_pt(i) = str2double(sub_pt{i}(3:end));
end

[C,IA,IB] = intersect(SUBNO_anim_pt, SUBNO);
TLESide_sel = TLESide(IB);

TLE_left = find(TLESide_sel == 'Left');

LI_anim_pt_val_left = LI_anim_pt_val(:,TLE_left,:);
LI_symb_pt_val_left = LI_symb_pt_val(:,TLE_left,:);

mLI_sub1 = squeeze(mean(LI_anim_pt_val_left,2));
mLI_sub2 = squeeze(mean(LI_symb_pt_val_left,2)); 

tag = 'anim vs. symb, pt (left)';

mLI_sub_pt_left = mLI_sub1 - mLI_sub2;

figure,
% subplot(3,1,2)
clear LI_val
for j=1:length(network_sel)
    hold on
    plot(mLI_sub_pt_left(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   500]);
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(tag)
ylabel('LI')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');

%%
mLI_sub_diff = mLI_sub_hc - mLI_sub_pt_left; tag = [mlabel,'; hc - pt-left'];

figure,
% subplot(3,1,3)
clear LI_val
for j=1:length(network_sel)
    hold on
    plot(mLI_sub_diff(network_sel(j),:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1100   500]);
end
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
title(tag)
ylabel('LI')
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');

set(gcf, 'Position', [1000   400   1000   900]);



