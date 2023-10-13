
%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 08/09/2023

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
clc
cfg = []; 
cfg.protocol = protocol;
cfg.datatag = 'wDICS_contrast_18_4';
cfg.datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';
cfg.BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_all_subjects';
S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);

%% Subject demog details
clc
cfg = []; cfg.subjs = S_data.subjs; 
cfg.sFiles = S_data.sFiles_32; 
sub_demog_data = ecpfunc_read_sub_demog_contrast(cfg);

%% TLE side (PT only)
patn_neuropsych_data = ecpfunc_read_patn_neuropsych();

%%
clc
cfg = [];
cfg.sub_demog_data = sub_demog_data;
cfg.patn_neuropsych_data = patn_neuropsych_data;
sub_TLE_sub_data = ecpfunc_read_sub_TLE_sub_contrast(cfg);

%% Inter-subject (group) averaging,
disp('1: Ctrl')
disp('2: Patn')

cfg = [];
cfg.sub_demog_data = sub_demog_data;
cfg.select_data = 1; S_data_hc = ecpfunc_select_contrast_data(cfg);
cfg.select_data = 2; S_data_pt = ecpfunc_select_contrast_data(cfg);

%%
cfg = []; cfg.strt = 0; cfg.spt = 2; cfg.overlap = 0.01; cfg.linterval = 0.3;
wi  = do_time_intervals(cfg);

%%
clc
glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat';
cfg = []; 
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
cfg.glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';
cfg.glass_atlas = glass_atlas;
Data_hcp_atlas = ecpfunc_hcp_atlas2(cfg);

net_sel_mutiple_label = Data_hcp_atlas.groups_labels';

%%
network_sel = [1:3,6:11];
colr = distinguishable_colors(length(network_sel));

%%
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/';
cd(data_save_dir)

%%
disp('1: threshold')
disp('2: counting')
disp('3: bootstrapping')
LI_method = input('LI_method sel:');
switch LI_method
    case 1
        mlabel = 'wDICS_threshold';
    case 2
        mlabel = 'wDICS_counting';
    case 3
        mlabel = 'wDICS_bsl_bootstrapping';
end

%%
outdir = fullfile(data_save_dir, mlabel, 'figs');
checkOrCreateDir(outdir);

cd(fullfile(data_save_dir, mlabel))

LI_hc = load('LI_Ctrl');
LI_pt = load('LI_Patn');

%% Power analysis
switch LI_method
    case 'threshold'
        run_power_analysis
end

%% Anim HC
mLI_sub1 = squeeze(nanmean(LI_hc.LI_sub,2)); tag = [mlabel, '; anim hc'];

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

plotData(mLI_sub_hc(network_sel, :),net_sel_mutiple_label(network_sel), colr, [mlabel, '; anim hc']);

%%
% - export figs
cfg = [];
cfg.outdir = outdir;
filename = tag;
cfg.filename = filename;
cfg.type = 'fig';
do_export_fig(cfg)

%% anim vs. symb - HC
% close all
mLI_sub1 = squeeze(nanmean(LI_hc.LI_sub,2));

tag = [mlabel,'; anim vs. symb, hc'];

clc
mLI_sub_hc = mLI_sub1;

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


% - export figs
cfg = [];
cfg.outdir = outdir;
filename = tag;
cfg.filename = filename;
cfg.type = 'fig';
do_export_fig(cfg)

%%
LI_pt_ID = S_data_pt.sFiles_subid;

mLI_sub1 = squeeze(mean(LI_pt.LI_sub,2));

tag = 'anim vs. symb, pt';

mLI_sub_pt = mLI_sub1;

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

% - export figs
cfg = [];
cfg.outdir = outdir;
filename = tag;
cfg.filename = filename;
cfg.type = 'fig';
do_export_fig(cfg)

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

% - export figs
cfg = [];
cfg.outdir = outdir;
cfg.filename = tag;
cfg.type = 'fig';
do_export_fig(cfg)

%%
patn_neuropsych_tle = ecpfunc_read_patn_neuropsych_tle();
TLESide = patn_neuropsych_tle.TLESide; SUBNO = patn_neuropsych_tle.SUBNO;

for i=1:length(LI_pt_ID)
    SUBNO_pt(i) = str2double(LI_pt_ID{i}(3:end));
end

[C,IA,IB] = intersect(SUBNO_pt, SUBNO);
TLESide_sel = TLESide(IB);

TLE_left = find(TLESide_sel == 'Left');

LI_pt_val_left = LI_pt.LI_sub(:,TLE_left,:);

mLI_sub1 = squeeze(mean(LI_pt_val_left,2));

tag = 'anim vs. symb, pt (left)';

mLI_sub_pt_left = mLI_sub1;

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

% - export figs
cfg = [];
cfg.outdir = outdir;
cfg.filename = tag;
cfg.type = 'fig';
do_export_fig(cfg)

%%
mLI_sub_diff = mLI_sub_hc - mLI_sub_pt_left; tag = [mlabel,'; hc - pt-left'];

figure,
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

% - export figs
cfg = []; cfg.outdir = outdir; cfg.filename = tag;
cfg.type = 'fig'; do_export_fig(cfg)

%% MEG vs. fMRI lat analysis (HC)
[sub_MF_hc,IA,IB] = intersect(S_data_hc.sFiles_subid, fmri_LIs.ID.language_Lateral);

mLI_sub1 = mean(squeeze(nanmean(LI_hc.LI_sub([1,2,6],:,:),1)),2);

tag = [mlabel,'; anim vs. symb, hc'];

% missing from fMRI 
difference = setdiff(S_data_hc.sFiles_subid, sub_MF_hc');
disp(difference');

%% MEG vs. fMRI lat analysis (PT)
clc
[sub_MF_pt,IA,IB] = intersect(LI_pt_ID, fmri_LIs.ID.language_Lateral');

LI_pt_val = LI_pt.LI_sub;

LI_pt_val_new = LI_pt_val(:,IA,:);
fmri_LIs_val = (fmri_LIs.val.language_Lateral(IB)); % Lateral regions was used.

% missing from fMRI 
difference = setdiff(LI_pt_ID, sub_MF_pt');
disp('missing from fMRI')
disp(difference');

%% MEG LI vs fMRI LI (language_Lateral)
% close all
clc

% Corr, MEG-fMRI
cfg = []; cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.ternary = 0;
cfg.thre = .2;
cfg.savefig = 1;
cfg.bf = 10;
cfg.outdir = outdir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new; 
cfg.fmri_LIs_val = fmri_LIs_val; 
% cfg.net_sel = [1,2,6];
cfg.net_sel = [11];
% cfg.net_sel = [6];
% cfg.net_sel = [1];
% cfg.net_sel = [2];
% [megLI_sub_pt, fmri_LIs_val, crr] = do_MEG_fMRI_corr(cfg);
[megLI_sub_pt, fmri_LIs_val, crr] = do_MEG_fMRI_corr_contrast(cfg);


%% MEG LI vs fMRI LI (Ternary language_Lateral)
close all
clc

cfg = [];
cfg.thre = .3; cfg.LI = fmri_LIs_val;
fmri_LIs_trn = do_ternary_classification(cfg);
size(fmri_LIs_trn);

% concordance, MEG-fMRI
cfg = [];
cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.ternary = 1;
cfg.savefig = 1;
cfg.outdir = outdir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new; 
cfg.fmri_LIs_val = fmri_LIs_trn; 
% cfg.net_sel = [1,2,6];
cfg.net_sel = [11];
cfg.thre = 0.15;
cfg.buffervalue = 4;
[megLIs_trn, fmri_LIs_trn] = do_MEG_fMRI_concordance_contrast(cfg);

% disp([megLIs_trn, fmri_LIs_trn])

idx = find(abs(megLIs_trn - fmri_LIs_trn) > 0);
sub_MF_pt(idx)

% length(megLIs_trn)

% (length(megLIs_trn) - length(idx))/length(megLIs_trn)

% figure, imagesc([megLIs_trn - fmri_LIs_trn])
% colorbar

%%
% close all
% clc
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other/')

cfg = [];
cfg.DataArray = [megLIs_trn, fmri_LIs_trn]; 
cfg.savefig = 1;
cfg.outdir = outdir; 
cfg.title = 'SD task, 58 PTs, ternary class';
do_plot_LIs(cfg)

cfg = [];
cfg.DataArray = [megLI_sub_pt, fmri_LIs_val]; 
cfg.savefig = 1;
cfg.outdir = outdir; 
cfg.title = 'SD task, 58 PTs, raw LIs';
do_plot_LIs(cfg)

size(fmri_LIs_val)

%%
[C, order] = confusionmat(megLIs_trn, fmri_LIs_trn);
% Visualize the confusion matrix
figure;
h = heatmap(C);
h.XDisplayLabels = {'-1','0', '1'};
h.YDisplayLabels = {'-1','0', '1'};
xlabel('MEG');
ylabel('fMRI');
title('Confusion Matrix');

% - export figs
cfg = []; cfg.outdir = outdir; cfg.filename = ['Confusion Matrix']; 
cfg.type = 'fig'; do_export_fig(cfg)

cd(outdir)

%% Corr. ROIs
% clc, close all

fmri_LIs_val = (fmri_LIs.val.language_Frontal(IB)); net_sel = 2;
% fmri_LIs_val = (fmri_LIs.val.language_Temporal(IB)); net_sel = 6;
fmri_LIs_val = (fmri_LIs.val.language_Lateral(IB)); net_sel = 11;


cfg = []; cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.thre = 0.1;
cfg.bf = 10;
cfg.ternary = 0;
cfg.savefig = 0;
cfg.outdir = outdir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new; 
cfg.fmri_LIs_val = fmri_LIs_val; cfg.net_sel = net_sel;
crr = do_MEG_fMRI_corr_contrast(cfg);

%% ROIs
clc
cfg = []; cfg.wi = wi;
cfg.ID = sub_MF_pt;
cfg.thre = 0.1;
cfg.bf = 10;
cfg.ternary = 0;
cfg.savefig = 0;
cfg.outdir = outdir;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new; 
cfg.fmri_LIs_val = fmri_LIs;
cfg.idx = IB;
% cfg.net_sel = net_sel;
crr = do_MEG_fMRI_corr_contrast_rois(cfg);

%% Corr.(tern)
clc, close all

fmri_LIs_val = (fmri_LIs.val.language_Frontal(IB)); net_sel = 2;
fmri_LIs_val = (fmri_LIs.val.language_Angular(IB)); net_sel = 1;
fmri_LIs_val = (fmri_LIs.val.language_Temporal(IB)); net_sel = 6;
fmri_LIs_val = (fmri_LIs.val.language_Lateral(IB)); net_sel = 11;

cfg = [];
cfg.thre = .3; cfg.LI = fmri_LIs_val;
fmri_LIs_trn = do_ternary_classification(cfg);
size(fmri_LIs_trn);

%- concordance (similarity)
cfg = [];
cfg.wi = wi;
cfg.thre = 0.15;
cfg.ternary = 1;
cfg.savefig = 0;
cfg.outdir = outdir;
cfg.ID = sub_MF_pt;
cfg.net_sel_mutiple_label = net_sel_mutiple_label;
cfg.LI_val = LI_pt_val_new;
% cfg.LI_symb_val = LI_symb_pt_val_new;
cfg.fmri_LIs_val = fmri_LIs_trn;
cfg.net_sel = net_sel;
cfg.buffervalue = 10;
conc = do_MEG_fMRI_concordance_contrast(cfg);

%%