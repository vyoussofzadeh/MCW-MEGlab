
%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 06/09/2023

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
% net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'BTLA'; 'VWFA'; 'Ang-Fron-Temp'};
net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Occipital'; 'Other'; 'PCingPrecun';'Temporal'; 'BTLA'; 'VWFA'};

%%
network_sel = [1:3,5:8];
colr = distinguishable_colors(length(network_sel));

%%
data_save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/LI_subs/group_8net_300ms';
cd(data_save_dir)

%%
cd(data_save_dir)

LI_anim_hc = load('LI_anim-hc');
LI_anim_pt = load('LI_anim-pt');

LI_symb_hc = load('LI_symb-hc');
LI_symb_pt = load('LI_symb-pt');

%% UPDATE POWER
cfg = []; cfg.d_in = LI_anim_hc; LI_anim_hc = do_readpow(cfg);
cfg = []; cfg.d_in = LI_symb_hc; LI_symb_hc = do_readpow(cfg);
cfg = []; cfg.d_in = LI_anim_pt; LI_anim_pt = do_readpow(cfg);
cfg = []; cfg.d_in = LI_symb_pt; LI_symb_pt = do_readpow(cfg);

LI = [];
LI.LI_anim_hc = LI_anim_hc;
LI.LI_symb_hc = LI_symb_hc;
LI.LI_anim_pt = LI_anim_pt;
LI.LI_symb_pt = LI_symb_pt;

[sub_pt,IA,IB] = intersect(S_data_anim_pt.sFiles_subid, S_data_symb_pt.sFiles_subid);

LI.LI_anim_pt.LI_sub = LI.LI_anim_pt.LI_sub(:,IA,:);
LI.LI_anim_pt.m_LI_max_sub = LI.LI_anim_pt.m_LI_max_sub(IA);
LI.LI_anim_pt.pow_sub = LI.LI_anim_pt.pow_sub(:,IA);
LI.LI_anim_pt.pow_LH = LI.LI_anim_pt.pow_LH(:,IA,:);
LI.LI_anim_pt.pow_RH = LI.LI_anim_pt.pow_RH(:,IA,:);

S_data_anim_pt.sFiles_in = S_data_anim_pt.sFiles_in(IA);
S_data_anim_pt.sFiles_subid = S_data_anim_pt.sFiles_subid(IA);

%% mean LI
clc, close all

S_in = [];
S_in{1} = S_data_anim_hc;
S_in{2} = S_data_symb_hc;
S_in{3} = S_data_anim_pt;
S_in{4} = S_data_symb_pt;

LI_in = [];
LI_in{1} = LI.LI_anim_hc;
LI_in{2} = LI.LI_symb_hc;
LI_in{3} = LI.LI_anim_pt;
LI_in{4} = LI.LI_symb_pt;

cfg = []; cfg.network_sel = [network_sel, 9];
cfg.wi = wi; cfg.colr = distinguishable_colors(length(cfg.network_sel));
cfg.net_label = [net_sel_mutiple_label; 'Ang-Fron-Temp'];
cfg.outdir = data_save_dir;
for j=1:length(S_in)
    mLI = squeeze(mean(LI_in{j}.LI_sub,2));
    cfg.mLI_sub = [mLI; mean(mLI([1,2,6],:),1)];
    cfg.S_data = S_in{j}.s_tag;
    do_plot_tLI(cfg)
end

%% neuropsych_tle
clc, close all
patn_neuropsych_tle = ecpfunc_read_patn_neuropsych_tle();

cfg = [];
cfg.patn_neuropsych_tle = patn_neuropsych_tle;
cfg.sub_pt = sub_pt;
[tle_IA, tle_IB, var_metric, var_name, TLE_idx] = do_extract_tle_data(cfg);

% - neuropsych_tle
cfg = [];
cfg.IA = tle_IA;
cfg.IB = tle_IB;
cfg.var_metric = var_metric;
cfg.TLE_idx = TLE_idx;
cfg.d_in = S_data_anim_pt; S_data_anim_pt_new_tle = do_update_tle_data(cfg);
cfg.d_in = S_data_symb_pt; S_data_symb_pt_new_tle = do_update_tle_data(cfg);

% -LI (& Pow)-pt_tle
cfg = [];
cfg.tle_IA = tle_IA;
cfg.LI = LI;
cfg.TLE_idx = TLE_idx;
LI_updated = do_update_tle_li_pow(cfg);

%%
clc
spec = [];
spec.network_sel = [network_sel, 9];
spec.colr = distinguishable_colors(length(spec.network_sel));
spec.net_label = [net_sel_mutiple_label; 'Ang-Fron-Temp'];
spec.outdir = data_save_dir;
spec.wi = wi;

%% LI
clc, close all
cfg = [];
cfg.d1 = LI_updated.LI_anim_hc.LI_sub;
cfg.d2 = LI_updated.LI_symb_hc.LI_sub;
cfg.spec = spec; cfg.spec.s_tag = 'LI anim vs. symb, hc';
LI_hc = do_tLI_diff(cfg);

cfg = [];
cfg.d1 = LI_updated.LI_anim_pt.LI_sub_tle;
cfg.d2 = LI_updated.LI_symb_pt.LI_sub_tle;
cfg.spec = spec; cfg.spec.s_tag = 'LI anim vs. symb, tle';
LI_pt = do_tLI_diff(cfg);

cfg = [];
cfg.d1 = LI_updated.LI_anim_pt.LI_sub_tle_left;
cfg.d2 = LI_updated.LI_symb_pt.LI_sub_tle_left;
cfg.spec = spec; cfg.spec.s_tag = 'LI anim vs. symb, pt, tle-left';
LI_tle_left = do_tLI_diff(cfg);

cfg = [];
cfg.d1 = LI_updated.LI_anim_pt.LI_sub_tle_right;
cfg.d2 = LI_updated.LI_symb_pt.LI_sub_tle_right;
cfg.spec = spec; cfg.spec.s_tag = 'LI anim vs. symb, pt, tle-right';
LI_tle_right = do_tLI_diff(cfg);

cfg = [];
cfg.network_sel = spec.network_sel;
cfg.wi = spec.wi;
cfg.colr = spec.colr;
cfg.net_label = spec.net_label;
cfg.outdir = spec.outdir;
cfg.mLI_sub = LI_tle_left - LI_tle_right;
cfg.S_data = 'LI tle-left vs tle-right';
do_plot_tLI(cfg)

cfg = [];
cfg.network_sel = spec.network_sel;
cfg.wi = spec.wi;
cfg.colr = spec.colr;
cfg.net_label = spec.net_label;
cfg.outdir = spec.outdir;
cfg.mLI_sub = LI_hc - LI_tle_left;
cfg.S_data = 'LI HC vs tle-left';
do_plot_tLI(cfg)

%% Pow - Left
clc, close all
cfg = [];
cfg.d1 = LI_updated.LI_anim_hc.pow_LH;
cfg.d2 = LI_updated.LI_symb_hc.pow_LH;
cfg.spec = spec; cfg.spec.s_tag = 'Pow anim vs. symb, hc';
Pow_hc = do_tLI_diff(cfg);

cfg = [];
cfg.d1 = LI_updated.LI_anim_pt.pow_sub_tle_LH;
cfg.d2 = LI_updated.LI_symb_pt.pow_sub_tle_LH;
cfg.spec = spec; cfg.spec.s_tag = 'Pow anim vs. symb, tle';
Pow_pt = do_tLI_diff(cfg);

cfg = [];
cfg.d1 = LI_updated.LI_anim_pt.pow_sub_tle_LH_left;
cfg.d2 = LI_updated.LI_symb_pt.pow_sub_tle_LH_left;
cfg.spec = spec; cfg.spec.s_tag = 'Pow anim vs. symb, pt, tle-left';
Pow_tle_left = do_tLI_diff(cfg);

cfg = [];
cfg.d1 = LI_updated.LI_anim_pt.pow_sub_tle_LH_right;
cfg.d2 = LI_updated.LI_symb_pt.pow_sub_tle_LH_right;
cfg.spec = spec; cfg.spec.s_tag = 'Pow anim vs. symb, pt, tle-right';
Pow_tle_right = do_tLI_diff(cfg);

cfg = [];
cfg.network_sel = spec.network_sel;
cfg.wi = spec.wi;
cfg.colr = spec.colr;
cfg.net_label = spec.net_label;
cfg.outdir = spec.outdir;
cfg.mLI_sub = Pow_tle_left - Pow_tle_right;
cfg.S_data = 'Pow tle-left vs tle-right';
do_plot_tLI(cfg)

cfg = [];
cfg.network_sel = spec.network_sel;
cfg.wi = spec.wi;
cfg.colr = spec.colr;
cfg.net_label = spec.net_label;
cfg.outdir = spec.outdir;
cfg.mLI_sub = Pow_hc - Pow_tle_left;
cfg.S_data = 'Pow HC vs tle-left';
do_plot_tLI(cfg)

%%
clc, close all
spec = [];
spec.wi = wi;
spec.plotflag = 0;
spec.sid = S_data_anim_hc.sFiles_subid;
spec.BS_dir = BS_data_dir;
spec.src = src;
spec.outdir = data_save_dir;
spec.net_sel_mutiple_label = {'Angular'; 'Frontal'; 'Temporal'};

cfg = [];
cfg.LI_sub = LI_updated.LI_anim_hc.LI_sub([1,2,6],:,:) - LI_updated.LI_symb_hc.LI_sub([1,2,6],:,:);
cfg.s1 = S_data_anim_hc; cfg.s2 = S_data_symb_hc;
spec.sid = S_data_anim_hc.sFiles_subid;
cfg.spec = spec; cfg.spec.stag = 'LI_sub_hc';
[LI_sub_hc, wi_sub_max_hc, sval_sub_hc] = do_maxLI_sourcemap_diff(cfg);

cfg = [];
cfg.LI_sub = LI_updated.LI_anim_pt.LI_sub([1,2,6],:,:) - LI_updated.LI_symb_pt.LI_sub([1,2,6],:,:);
cfg.s1 = S_data_anim_pt; cfg.s2 = S_data_symb_pt;
spec.sid = S_data_anim_pt.sFiles_subid;
cfg.spec = spec; cfg.spec.stag = 'LI_sub_pt';
[LI_sub_pt, wi_sub_max_pt, sval_sub_pt] = do_maxLI_sourcemap_diff(cfg);

cfg = [];
cfg.LI_sub = LI_updated.LI_anim_pt.LI_sub_tle([1,2,6],:,:) - LI_updated.LI_symb_pt.LI_sub_tle([1,2,6],:,:);
cfg.s1 = S_data_anim_pt_new_tle.d_tle; cfg.s2 = S_data_symb_pt_new_tle.d_tle;
spec.sid = S_data_anim_pt_new_tle.d_tle.sFiles_subid;
cfg.spec = spec; cfg.spec.stag = 'LI_sub_tle';
[LI_sub_tle, wi_sub_max_tle, sval_sub_tle] = do_maxLI_sourcemap_diff(cfg);

cfg = [];
cfg.LI_sub = LI_updated.LI_anim_pt.LI_sub_tle_left([1,2,6],:,:) - LI_updated.LI_symb_pt.LI_sub_tle_left([1,2,6],:,:);
cfg.s1 = S_data_anim_pt_new_tle.d_left; cfg.s2 = S_data_symb_pt_new_tle.d_left;
spec.sid = S_data_anim_pt_new_tle.d_left.sFiles_subid;
cfg.spec = spec; cfg.spec.stag = 'LI_sub_ltle';
[LI_sub_ltle, wi_sub_max_ltle, sval_sub_ltle] = do_maxLI_sourcemap_diff(cfg);

cfg = [];
cfg.LI_sub = LI_updated.LI_anim_pt.LI_sub_tle_right([1,2,6],:,:) - LI_updated.LI_symb_pt.LI_sub_tle_right([1,2,6],:,:);
cfg.s1 = S_data_anim_pt_new_tle.d_right; cfg.s2 = S_data_symb_pt_new_tle.d_right;
spec.sid = S_data_anim_pt_new_tle.d_right.sFiles_subid;
cfg.spec = spec; cfg.spec.stag = 'LI_sub_rtle';
[LI_sub_rtle, wi_sub_max_rtle, sval_sub_rtle] = do_maxLI_sourcemap_diff(cfg);

%% Corr analysis
close all
cfg = [];
cfg.d_in = mean(LI_sub_ltle,2);
% cfg.d_in = max(LI_sub_ltle');
cfg.var_name = var_name; cfg.metric = S_data_anim_pt_new_tle.d_left.all_metrics;
do_corr_analysis(cfg);

cfg = [];
cfg.d_in = mean(LI_sub_rtle,2);
cfg.d_in = max(LI_sub_rtle');
cfg.var_name = var_name; cfg.metric = S_data_anim_pt_new_tle.d_right.all_metrics;
do_corr_analysis(cfg)

cfg = [];
cfg.d_in = mean(LI_sub_tle,2);
cfg.d_in = max(LI_sub_tle');
cfg.var_name = var_name; cfg.metric = S_data_anim_pt_new_tle.d_tle.all_metrics;
do_corr_analysis(cfg)

%%
spec = [];
spec.wi = wi;
spec.outdir = fullfile(outdir,'group_8net');
spec.net_label = net_sel_mutiple_label;
spec.net_sel = [1:2,6:8];

cfg = []; cfg.wi = spec.wi; cfg.savefig = 1;
cfg.outdir = spec.outdir;
cfg.net_sel_mutiple_label = spec.net_label ;
cfg.network_sel = spec.net_sel;

cfg.LI_sub = LI_updated.LI_anim_hc.LI_sub;
cfg.S_data_sel = S_data_anim_pt;
do_plot_group_lat(cfg);

cfg.LI_sub = LI_updated.LI_anim_pt.LI_sub;
cfg.S_data_sel = S_data_anim_pt;
do_plot_group_lat(cfg);

cfg.LI_sub = LI_updated.LI_anim_pt.LI_sub_tle_left;
cfg.S_data_sel = S_data_anim_pt_new_tle.d_left;
do_plot_group_lat(cfg);

cfg.LI_sub = LI_updated.LI_anim_pt.LI_sub_tle;
cfg.S_data_sel = S_data_anim_pt_new_tle.d_tle;
do_plot_group_lat(cfg);

%%
figure,
cfg = [];
cfg.view = [-180,-90;0,90;-90,0; 90,0];
cfg.color = viridis(256);
cfg.position = [800   800   1000   300];
cfg.title = '';
cfg.alpha = 1; cfg.coor = [];
cfg.surf = src;

tmp = m_source_hc; tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); m_source_hc_norm = tmp;
tmp = m_source_tle_left; tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); m_source_tle_left_norm = tmp;
tmp = m_source_tle_right; tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); m_source_tle_right_norm = tmp;

cfg.d_in = (m_source_hc_norm - m_source_tle_left_norm);
cfg.d_in = (m_source_hc_norm - m_source_tle_right_norm);
% cfg.d_in = (m_source_tle_left_norm - m_source_tle_right_norm);
do_surfplot(cfg);

%%
tmp = sval_sub_hc; tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); m_source_hc_norm = tmp;
tmp = sval_sub_ltle; tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); m_source_tle_left_norm = tmp;
tmp = sval_sub_rtle; tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); m_source_tle_right_norm = tmp;

figure,
cfg = [];
cfg.view = [-180,-90;0,90;-90,0; 90,0];
cfg.color = viridis(256);
cfg.position = [800   800   1000   300];
cfg.title = '';
cfg.alpha = 1; cfg.coor = [];
% m_source = sval_sub_rtle';
% tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:)));
cfg.surf = src;
cfg.d_in = (m_source_hc_norm - m_source_tle_left_norm);
cfg.d_in = (m_source_tle_left_norm);
% cfg.d_in = (m_source_hc_norm - m_source_tle_right_norm);% - sval_sub_rtle;
do_surfplot(cfg);

