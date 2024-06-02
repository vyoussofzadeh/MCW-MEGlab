% MEG Time-Frequency Analysis Pipeline
% Written by MCW Group, Vahab Youssofzadeh <vyoussofzadeh@mcw.edu>
% Latest Update: 05/11/2024

%% Clear workspace and figures
clc; close all

%% Setup path
ft_path = '/opt/matlab_toolboxes/ft_packages/fieldtrip_latest';
addpath(ft_path);
ft_defaults

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/psb_pilot/functions')

%% Load data
Data = load('/group/prishah/LanguageMEG/HC005/meg/baseline/all_sensdata_clean_pstm.mat');

%% Preprocessing, band-pass filtering
cfg = [];
cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.hpfiltord = 3;
cfg.hpfreq = 0.1;
cfg.lpfreq = 55;
cfg.channel = {'MEG'};
f_data = ft_preprocessing(cfg, Data.alldata);

%% Processing, trial selection
cfg = [];
cfg.trials = f_data.trialinfo == 4;
D_4 = ft_selectdata(cfg, Data.alldata);

cfg.trials = f_data.trialinfo == 8;
D_8 = ft_selectdata(cfg, Data.alldata);

%% Processing, interval selection
cfg = [];
cfg.latency = [-1, 4];
D_4_interval = ft_selectdata(cfg, D_4);
D_8_interval = ft_selectdata(cfg, D_8);

%% Preprocessing, notch filtering
cfg = [];
cfg.subj = '005';
nD_4_interval = do_notch(cfg, D_4_interval);
nD_8_interval = do_notch(cfg, D_8_interval);

%% Processing, freq analysis
cfg = []; 
cfg.savefile = []; 
cfg.saveflag = 0;
cfg.foilim = [1 45]; 
cfg.plotflag  = 1;
cfg.tapsmofrq = 4; 
cfg.taper = 'hanning';
do_fft(cfg, nD_4_interval);
do_fft(cfg, nD_8_interval);

%% Processing, layout selection
cfg = []; cfg.layout = 'neuromag306mag.lay'; lay = ft_prepare_layout(cfg);

%% Processing, time-freq analysis
cfg = []; cfg.layout = lay; cfg.subj = []; cfg.baseline = [-0.5,0];
cfg.fmax = 45; cfg.title = 'Time-Freq';
[t_max_4,f_max_4, tfr_4]  = do_tfr_analysis(cfg, nD_4_interval);
[t_max_8,f_max_8, tfr_8]  = do_tfr_analysis(cfg, nD_8_interval);

%% Processing, grand mean
aD_4 = do_ave(nD_4_interval);
aD_8 = do_ave(nD_8_interval);

cfg = []; cfg.savefile = []; cfg.saveflag = 2; cfg.lay  = lay; do_ave_plot(cfg, aD_4);
cfg = []; cfg.savefile = []; cfg.saveflag = 2; cfg.lay  = lay; do_ave_plot(cfg, aD_8);

