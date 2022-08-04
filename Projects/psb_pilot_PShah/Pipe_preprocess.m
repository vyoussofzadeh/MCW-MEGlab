%% The psb_pilot

% MEG (pre)-processing pipeline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 05/10/2022

clear; clc, close('all'); warning off

%% Paths
restoredefaultpath
script_path = '/data/MEG/Research/psb_pilot';
addpath(genpath(script_path));

%- Input dir
indir = '/data/MEG/Research/psb_pilot/ss_pilot/tsss';
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

%% Reading data, step 1: Layout & sensor location
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);

%
subj = 'pilot1';
datafile = fullfile(indir,'STM_Block1_raw.fif');
sens = ft_read_sens(datafile);
sens = ft_convert_units(sens,'mm');

%% Reading data, step 2: read event
% event = ft_read_event(datafile);
hdr = ft_read_header(datafile); %read header information
Fsample = hdr.Fs;

tt = linspace(1, length(detTrig)/Fsample, length(detTrig));

Index = strfind(hdr.label,{'STI101'});
Index = find(not(cellfun('isempty',Index)));

if isempty(Index)
    error('stimulus data, STI101, is missing from the data header, ');
end

detTrig=ft_read_data(datafile,'chanindx',Index,'header',hdr,'eventformat','neuromag_fif','dataformat','neuromag_fif');
detTrig = (detTrig - min(detTrig));
detTrig=bitand(detTrig,255);

figure,
subplot 211
plot(tt,detTrig)
title('detTrig')

detResp=ft_read_data(datafile,'chanindx',Index,'header',hdr,'eventformat','digital trigger','dataformat','neuromag_fif');
detResp = detResp - detTrig;
detResp = (detResp - min(detResp));
resp = (detResp - mean(detResp))./max(detTrig);

%     figure,
subplot 212
plot(tt,detResp)
title('detResp')

figure,
plot(tt,detTrig)
hold on
plot(tt,resp, 'r')
title('All triggers')

%% Preprocessing, step 1: Read & Band-pass filter
cfg = [];
cfg.datafile = datafile;
cfg.prestim = 1;
cfg.poststim = 3;
cfg.eventvalue = 1; % this can be changed based on the event type (see event var)
cfg.hpfreq  = 0.1;
cfg.lpfreq = 70;
cfg.dftfreq = 60; % or [60 120 180];
f_data = do_ft_preprocess(cfg);

%% Preprocessing, step 2: Rejecting bad data
cfg = [];
cfg.pflag = 1; % yes:1, No:2
cfg.saveflag = 0; % yes:1, No:2
cfg.savepath = [];
cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)]; 
cfg.rejectpercentage = .95;
cfg.method = 'auto'; % 'manual'
[r_data,report] = do_rejectdata(cfg, f_data);

%% Preprocessing, step 3: ICA
cfg = [];
cfg.lay = lay;
cfg.subj = subj;
cfg.n = 20;
cfg.allpath = allpath;
cfg.savefig = 1;
data_ica = do_ica(cfg, r_data);

%% Preprocessing, step 4: Notch filtering
cfg = [];
cfg.subj = subj;
data_clean = do_notch(cfg, data_ica);

%% Preprocessing, step 5: Save data
disp('saving data ...')
cd(outdir)
save(['dataclean_', subj], 'data_clean', '-v7.3');

%% Processing, step 1: Freq analysis
cfg = [];
cfg.savefile = [];
cfg.saveflag = 0;
cfg.foilim = [2 50];
cfg.plotflag  = 1;
cfg.tapsmofrq       = 4;
cfg.taper    = 'hanning';
do_fft(cfg, data_clean);

%% Processing, step 2: Time-Freq analysis
cfg = [];
cfg.layout = lay;
cfg.subj = subj;
cfg.baseline = [-0.3,0];
cfg.fmax = 50;
cfg.title = 'Time-Freq';
[t_max,f_max, tfr]  = do_tfr_analysis(cfg, data_clean);

%% Processing, step 3: Grand Mean
a_data = do_ave(data_clean);
cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.lay  = lay;
do_ave_plot(cfg, a_data);

%% Processing, step 4: Time-locked
% toi = [-0.3,0; 0,3]; ep_data = do_epoch(datain, toi);
