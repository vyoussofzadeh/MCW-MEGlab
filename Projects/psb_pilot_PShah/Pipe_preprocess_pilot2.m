%% The psb_pilot

% MEG (pre)-processing pipeline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 05/24/2022

clear; clc, close('all'); warning off

%% Flag
flag = [];
flag.plottriggers = 1;

%% Paths
restoredefaultpath
script_path = '/data/MEG/Research/psb_pilot';
addpath(genpath(script_path));

%- Input dir
indir = '/data/MEG/Research/psb_pilot/ss_pilot_2/220512/tsss';
%- Output dir
outdir = '/data/MEG/Research/psb_pilot/FT';

%
ft_path = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_20190419');
addpath(ft_path);
ft_defaults

allpath = [];
allpath.ft_path18 = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_041718'); % needed for IC plotting
allpath.ft_path = ft_path;
allpath.script_path = script_path;

%% Reading data, Layout & sensor location
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);

%
subj = 'pilot2';
datafile = fullfile(indir,'STM_Block1_raw.fif');
sens = ft_read_sens(datafile);
sens = ft_convert_units(sens,'mm');

%% Reading trigger info
cfg = [];
cfg.plot = 1;
[detResp,detTrig] = do_readtriggers(cfg, datafile);
    
%- THIS IS OPTIOPNAL
evt = ft_read_event(datafile);

%% Preprocessing, Read & Band-pass filter
cfg = [];
cfg.datafile = datafile;
cfg.prestim = 5;%1; in seconds
cfg.poststim = 0;%3; in seconds
cfg.eventvalue = [4 8];%[1, 2, 4, 8]; % this can be changed based on the event type (see event var)
cfg.hpfreq  = 0.1;
cfg.lpfreq = 70;
cfg.dftfreq = 60; % or [60 120 180];
f_data = do_ft_preprocess(cfg);

%% Preprocessing, Rejecting bad data
cfg = [];
cfg.pflag = 1; % yes:1, No:2
cfg.saveflag = 0; % yes:1, No:2
cfg.savepath = [];
cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)];
cfg.rejectpercentage = .95;
cfg.method = 'auto'; % 'manual'
[r_data,report] = do_rejectdata(cfg, f_data);

%% Preprocessing, ICA
cfg = [];
cfg.lay = lay;
cfg.subj = subj;
cfg.n = 20;
cfg.allpath = allpath;
cfg.savefig = 1;
data_ica = do_ica(cfg, r_data);

%% Preprocessing, Notch filtering
cfg = [];
cfg.subj = subj;
data_clean = do_notch(cfg, data_ica);

%% Preprocessing, Save data
disp('saving data ...')
cd(outdir)
save(['dataclean_', subj], 'data_clean', '-v7.3');

%% Processing, Freq analysis
cfg = [];
cfg.savefile = [];
cfg.saveflag = 0;
cfg.foilim = [2 50];
cfg.plotflag  = 1;
cfg.tapsmofrq       = 4;
cfg.taper    = 'hanning';
do_fft(cfg, data_clean);

%% Processing, Time-Freq analysis
cfg = [];
cfg.layout = lay;
cfg.subj = subj;
cfg.baseline = [-0.3,0];
cfg.fmax = 50;
cfg.title = 'Time-Freq';
[t_max,f_max, tfr]  = do_tfr_analysis(cfg, data_clean);

%% Processing, Grand Mean
a_data = do_ave(data_clean);
cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.lay  = lay;
do_ave_plot(cfg, a_data);

%% Processing, Time-locked
% toi = [-0.3,0; 0,3]; ep_data = do_epoch(datain, toi);

%% DATA SELECT
idx.cond1 = find(data_clean.trialinfo == 1);
idx.cond2 = find(data_clean.trialinfo == 2);
idx.cond3 = find(data_clean.trialinfo == 4);
idx.cond4 = find(data_clean.trialinfo == 8);

cfg = []; cfg.trials = idx.cond1; data_cond1 = ft_selectdata(cfg, data_clean);
cfg = []; cfg.trials = idx.cond2; data_cond2 = ft_selectdata(cfg, data_clean);
cfg = []; cfg.trials = idx.cond3; data_cond3 = ft_selectdata(cfg, data_clean);
cfg = []; cfg.trials = idx.cond4; data_cond4 = ft_selectdata(cfg, data_clean);

