%% The psb_pilot

% MEG (pre)-processing pipeline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 05/25/2022 by Shah-Basak, Priyanka

clear; clc, close('all'); warning off

%% Flag
flag = [];
flag.plottriggers = 1;

%% Paths
restoredefaultpath
maindir = '/group/prishah/LanguageMEG/psb_pilot';
ssid = 'ss_pilot_2';
datcol = '220512';

script_path = maindir;%'/data/MEG/Research/psb_pilot';
addpath(genpath(script_path));
%- Input dir
indir = fullfile(maindir, ssid, datcol, 'tsss');%'/data/MEG/Research/psb_pilot/ss_pilot_2/220512/tsss';
%- Output dir
outdir = fullfile(maindir, 'FT');

%
ft_path = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_20190419');
addpath(ft_path);
ft_defaults

allpath = [];
allpath.ft_path18 = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_041718'); % needed for IC plotting
allpath.ft_path = ft_path;
allpath.script_path = script_path;

%% Reading Layout & sensor location
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);
subj = 'pilot2';

%% Reading data, Layout & sensor location, Preprocessing, saving data
data_clean = [];
for bb = 1:6
    close all;
    datafile = fullfile(indir,['STM_Block', num2str(bb), '_raw.fif']);
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');

    % Reading trigger info
    cfg = [];
    cfg.plot = 1;
    [detResp,detTrig] = do_readtriggers(cfg, datafile);

    %- THIS IS OPTIOPNAL
    evt = ft_read_event(datafile);

    % Preprocessing, Read & Band-pass filter
    cfg = [];
    cfg.datafile = datafile;
    cfg.prestim = 5;%1; in seconds
    cfg.poststim = 0;%3; in seconds
    cfg.eventvalue = [4 8];%[1, 2, 4, 8]; % this can be changed based on the event type (see event var)
    cfg.hpfreq  = 0.1;
    cfg.lpfreq = 70;
    cfg.dftfreq = 60; % or [60 120 180];
    f_data = do_ft_preprocess(cfg);

    % Preprocessing, Rejecting bad data
    cfg = [];
    cfg.pflag = 1; % yes:1, No:2
    cfg.saveflag = 0; % yes:1, No:2
    cfg.savepath = [];
    cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)];
    cfg.rejectpercentage = .95;
    cfg.method = 'auto'; % 'manual'
    [r_data,report] = do_rejectdata(cfg, f_data);

    % Preprocessing, ICA
    cfg = [];
    cfg.lay = lay;
    cfg.subj = subj;
    cfg.n = 20;
    cfg.allpath = allpath;
    cfg.savefig = 1;
    data_ica = do_ica(cfg, r_data);

    % Preprocessing, Notch filtering
    %note psb padding changed to 5 from 4
    cfg = [];
    cfg.subj = subj;
    data_clean{bb} = do_notch(cfg, data_ica);
end
%save data
disp('saving data ...')
if ~exist(outdir,'dir')
    mkdir(outdir)
end
cd(outdir)
save(['dataclean_pstm_', subj], 'data_clean', '-v7.3');
%
%% Combine all blocks
%load(['data_clean_', subj])
alldata = data_clean{1};
fieldnames(alldata)

alldata.sampleinfo = cat(1, data_clean{1}.sampleinfo, ...
    data_clean{2}.sampleinfo, ...
    data_clean{3}.sampleinfo, ...
    data_clean{4}.sampleinfo, ...
    data_clean{5}.sampleinfo, ...
    data_clean{6}.sampleinfo);

alldata.trialinfo = cat(1, data_clean{1}.trialinfo, ...
    data_clean{2}.trialinfo, ...
    data_clean{3}.trialinfo, ...
    data_clean{4}.trialinfo, ...
    data_clean{5}.trialinfo, ...
    data_clean{6}.trialinfo);

alldata.trial = cat(2, data_clean{1}.trial, ...
    data_clean{2}.trial, ...
    data_clean{3}.trial, ...
    data_clean{4}.trial, ...
    data_clean{5}.trial, ...
    data_clean{6}.trial);

alldata.time = cat(2, data_clean{1}.time, ...
    data_clean{2}.time, ...
    data_clean{3}.time, ...
    data_clean{4}.time, ...
    data_clean{5}.time, ...
    data_clean{6}.time)
save(['all_dataclean_pstm_', subj], 'alldata', '-v7.3');
%% Processing, Freq analysis
cfg = [];
cfg.savefile = [];
cfg.saveflag = 0;
cfg.foilim = [2 50];
cfg.plotflag  = 1;
cfg.tapsmofrq = 4;
cfg.taper    = 'hanning';
do_fft(cfg, alldata);

%% Processing, Time-Freq analysis
cfg = [];
cfg.layout = lay;
cfg.subj = subj;
cfg.baseline = [-5,-4.7];
cfg.fmax = 50;
cfg.title = 'Time-Freq';
[t_max,f_max, tfr]  = do_tfr_analysis(cfg, alldata);

%% Processing, Grand Mean
a_data = do_ave(alldata);
cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.lay  = lay;
do_ave_plot(cfg, a_data);

%% Processing, Time-locked
% toi = [-0.3,0; 0,3]; ep_data = do_epoch(datain, toi);

%% DATA SELECT
%idx.cond1 = find(data_clean.trialinfo == 1);
%idx.cond2 = find(data_clean.trialinfo == 2);
idx.cond3 = find(alldata.trialinfo == 4);
idx.cond4 = find(alldata.trialinfo == 8);

%cfg = []; cfg.trials = idx.cond1; data_cond1 = ft_selectdata(cfg, data_clean);
%cfg = []; cfg.trials = idx.cond2; data_cond2 = ft_selectdata(cfg, data_clean);
cfg = []; cfg.trials = idx.cond3; data_cond3 = ft_selectdata(cfg, alldata);
cfg = []; cfg.trials = idx.cond4; data_cond4 = ft_selectdata(cfg, alldata);

