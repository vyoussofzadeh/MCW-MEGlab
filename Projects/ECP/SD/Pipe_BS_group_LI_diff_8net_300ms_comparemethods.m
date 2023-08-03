
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

%% Read fMRI lats
fmri_LIs = ecpfunc_read_fmri_lat();

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
net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; ...
    'PCingPrecun';'Temporal'; 'BTLA'; 'VWFA'};

cfg = []; Data_hcp_atlas = ecpfunc_hcp_atlas2(cfg);
net_sel_mutiple_label = Data_hcp_atlas.groups_labels';

%%
network_sel = [1:3,6:10];
colr = distinguishable_colors(length(network_sel));

%%
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/LI_subs/group_8net_300ms';
cd(data_save_dir)

%%
disp('1: threshold')
disp('2: counting')
disp('3: bootstrapping')
LI_method = input('LI_method sel:');
switch LI_method
    case 1
        mlabel = 'threshold_1';
    case 2
        mlabel = 'counting';
    case 3
        mlabel = 'bootstrapping';
end
%%
cd(fullfile(data_save_dir, mlabel))

LI_anim_hc = load('LI_anim-hc');
LI_anim_pt = load('LI_anim-pt');

LI_symb_hc = load('LI_symb-hc');
LI_symb_pt = load('LI_symb-pt');

%% Power analysis
close all

switch LI_method
    case 'threshold'
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
        
end

%% Anim HC
% close all
mLI_sub1 = squeeze(nanmean(LI_anim_hc.LI_sub,2)); tag = [mlabel, '; anim hc'];

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
mLI_sub1 = squeeze(nanmean(LI_anim_hc.LI_sub,2));
mLI_sub2 = squeeze(nanmean(LI_symb_hc.LI_sub,2)); tag = [mlabel,'; anim vs. symb, hc'];

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

LI_pt_ID = S_data_anim_pt.sFiles_subid(IA);

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

%% MEG vs. fMRI lat analysis (HC)
[sub_MF_hc,IA,IB] = intersect(S_data_anim_hc.sFiles_subid, fmri_LIs.ID.language_Lateral);

mLI_sub1 = mean(squeeze(nanmean(LI_anim_hc.LI_sub([1,2,6],:,:),1)),2);
mLI_sub2 = mean(squeeze(nanmean(LI_symb_hc.LI_sub([1,2,6],:,:),1)),2);
mLI_sub_hc = mLI_sub1 - mLI_sub2;

tag = [mlabel,'; anim vs. symb, hc'];

figure, plot(mLI_sub_hc(IA), str2double(fmri_LIs.val.language_Lateral(IB)),'*')

%%
[sub_MF_pt,IA,IB] = intersect(LI_pt_ID, fmri_LIs.ID.language_Lateral);

%% MEG LI vs fMRI LI (language_Lateral)
close all
clc

LI_anim_pt_val_new = LI_anim_pt_val(:,IA,:);
LI_symb_pt_val_new = LI_symb_pt_val(:,IA,:);

fmri_LIs_val = str2double(fmri_LIs.val.language_Lateral(IB));
% fmri_LIs_val = str2double(fmri_LIs.val.semantic_Lateral(IB));


% san_check = [LI_pt_ID(IA); fmri_LIs.ID.language_Lateral(IB)']';

cr = [];
for i=1:length(wi)
    %     for j=1:size(LI_anim_pt_val_new,2)
    mLI_sub1 = mean(LI_anim_pt_val_new([1,2,6],:,i));
    mLI_sub2 = mean(LI_symb_pt_val_new([1,2,6],:,i));
    mLI_sub_pt = (mLI_sub1 - mLI_sub2)';
%     mLI_sub_pt = mLI_sub1
    cr(i,:) = corr2(mLI_sub_pt, fmri_LIs_val);
    %     end
end

figure, plot(mean(wi'),cr), title('Lateral');
ylabel('LIs corr (MEG , fMRI)')

mLI_sub1 = mean(LI_anim_pt_val_new([1,2,6],:,idx));
mLI_sub2 = mean(LI_symb_pt_val_new([1,2,6],:,idx));
mLI_sub_pt = (mLI_sub1 - mLI_sub2)';

figure, bar([mLI_sub_pt, fmri_LIs_val])
figure, plot(mLI_sub_pt, fmri_LIs_val,'*')

corr2(mLI_sub_pt, fmri_LIs_val)

%%
close all
y = cr;
t =  mean(wi');
% Find the global peak and its index
[max_val, idx] = max(y);

% Define an interval of 10 samples centered around the peak
interval_length = 15;
interval_start = max(1, idx - interval_length/2);
interval_end = min(length(y), idx + interval_length/2);

% Get the interval containing the peak
peak_interval = y(interval_start:interval_end);

% Plot the original signal
plot(t, y)
hold on

% Highlight the interval containing the peak
plot(t(interval_start:interval_end), peak_interval, 'r')
hold off

%% Ternary classification
close all

%- MEG
% cfg = []; cfg.thre = 10;
% cfg.LI = LI_anim_pt_val_new; LI_anim_pt_trn = do_ternary_classification(cfg);
% cfg.LI = LI_symb_pt_val_new; LI_symb_pt_trn = do_ternary_classification(cfg);
% 
% size(LI_symb_pt_trn);

% figure, plot(squeeze(LI_anim_pt_val_new(1,1,:)))
% figure, plot(squeeze(LI_anim_pt_trn(1,1,:)),'r')

%- fMRI
cfg = [];
cfg.thre = 30;
cfg.LI = fmri_LIs_val;
fmri_LIs_trn = do_ternary_classification(cfg);

size(fmri_LIs_trn);

% figure, bar((fmri_LIs_val))
% figure, bar(fmri_LIs_trn,'r')

%%
%- MEG vs. fMRI (ternary)
dd = [];
for i=1:length(wi)
    %     for j=1:size(LI_anim_pt_val_new,2)
    mLI_sub1 = mean(LI_anim_pt_val_new([1,2,6],:,i));
    mLI_sub2 = mean(LI_symb_pt_val_new([1,2,6],:,i));
    %     mLI_sub_pt = round(mLI_sub1 - mLI_sub2)';
    mLI_sub_pt = (mLI_sub1 - mLI_sub2)';
    
    cfg = []; cfg.thre = 10;
    cfg.LI = mLI_sub_pt; mLI_sub_pt_trn = do_ternary_classification(cfg);
    
    %     mLI_sub_pt = mLI_sub1';
    dd(i,:) = (mLI_sub_pt_trn .* fmri_LIs_trn);
    %     end
end
figure, plot(mean(wi'),mean(dd,2)), title('Lateral');
ylabel('LIs (MEG * fMRI)')

%% Corr.
clc, close all

LI_anim_pt_val_new = LI_anim_pt_val(:,IA,:);
LI_symb_pt_val_new = LI_symb_pt_val(:,IA,:);

fmri_LIs_val = str2double(fmri_LIs.val.language_Frontal(IB)); net_sel = 2;
fmri_LIs_val = str2double(fmri_LIs.val.language_Angular(IB)); net_sel = 1;
% fmri_LIs_val = str2double(fmri_LIs.val.language_Temporal(IB)); net_sel = 6;
% 

cfg = []; cfg.wi = wi;
cfg.LI_anim_val = LI_anim_pt_val_new; cfg.LI_symb_val = LI_symb_pt_val_new;
cfg.fmri_LIs_val = fmri_LIs_val; cfg.net_sel = net_sel;
crr = do_MEG_fMRI_corr(cfg);


figure, plot(mean(wi'),crr), title([net_sel_mutiple_label{net_sel}]);
ylabel('LIs corr (MEG , fMRI)')

[mx, idx] = max(crr);

mLI_sub1 = (LI_anim_pt_val_new(net_sel,:,idx));
mLI_sub2 = (LI_symb_pt_val_new(net_sel,:,idx));
mLI_sub_pt = (mLI_sub1 - mLI_sub2)';

figure, bar([mLI_sub_pt, fmri_LIs_val])
figure, plot(mLI_sub_pt, fmri_LIs_val,'*')

corr2(mLI_sub_pt, fmri_LIs_val)

%% Tern, concordance (similarity)
clc, close all

LI_anim_pt_val_new = LI_anim_pt_val(:,IA,:);
LI_symb_pt_val_new = LI_symb_pt_val(:,IA,:);

% fmri_LIs_val = str2double(fmri_LIs.val.language_Frontal(IB)); net_sel = 2;
fmri_LIs_val = str2double(fmri_LIs.val.language_Angular(IB)); net_sel = 1;
% fmri_LIs_val = str2double(fmri_LIs.val.language_Temporal(IB)); net_sel = 6;

cfg = [];
cfg.wi = wi;
cfg.LI_anim_val = LI_anim_pt_val_new;
cfg.LI_symb_val = LI_symb_pt_val_new;
cfg.fmri_LIs_val = fmri_LIs_val;
cfg.net_sel = net_sel;
conc = do_MEG_fMRI_concordance(cfg);

figure, plot(mean(wi'),mean(conc,2)), title([net_sel_mutiple_label{net_sel}]);
ylabel('LIs corr (MEG , fMRI)')














