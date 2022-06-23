%% The bgros MEG pilot data

% MEG analysis pipeline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 06/09/2022

clear; clc, close('all'); warning off

%% Analysis flags
flag = [];
flag.plottriggers = 1;
flag.anatomy_check = 1;
flag.dicsanalysis = 0;
flag.connanalysis = 1;

%% Paths
addpath('./run')
Run_setpath

%% Sub specifications
subj = 'pilot1';
disp(subj)
run = input('enter run of intersts, e.g, 1: ');
switch run
    case {6, 7}
        run = '6_7';
end
condval = input('enter condition of intersts, e.g, 1, 2, 4, 8: ');
datafilename = ['Task_run', num2str(run), '_raw.fif'];
datafile = fullfile(indir, datafilename);
subdir = fullfile(outdir, subj,'preprocess'); % output dir
load(fullfile(subdir,['data_run', num2str(run), '_cond',num2str(condval), '.mat']))

savedir = fullfile(outdir, subj,'process',['run',num2str(run)]); % output dir
if exist(savedir, 'file') == 0, mkdir(savedir); end
cd(savedir)

%% Processing, Freq analysis
cfg = []; cfg.savefile = []; cfg.saveflag = 0;
cfg.foilim = [2 50]; cfg.plotflag  = 1;
cfg.tapsmofrq = 4; cfg.taper = 'hanning';
do_fft(cfg, data_clean);

%% Processing, Time-Freq analysis
cfg = []; cfg.layout = 'neuromag306mag.lay'; lay = ft_prepare_layout(cfg);

cfg = []; cfg.layout = lay; cfg.subj = subj; cfg.baseline = [-0.3,0];
cfg.fmax = 50; cfg.title = 'Time-Freq';
[t_max,f_max, tfr]  = do_tfr_analysis(cfg, data_clean);

%% Processing, Grand mean
a_data = do_ave(data_clean);

cfg = []; cfg.savefile = []; cfg.saveflag = 2; cfg.lay  = lay; do_ave_plot(cfg, a_data);

%% Processing, time-locked analysis
disp('enter time interval of interests, e.g, [-0.3,0; 0.4,1]')
toi  = input(':');
ep_data = do_epoch(data_clean, toi);

%% Anatomy
mridir = '/group/bgross/work/CIDMEG/analysis/anatomy/Sub001/T1';
mrifile = 'sub-CIDFMRI001_ses-1_acq-mprage_T1w.nii.gz';
Run_anatomy

%% Source analysis, time-domain (LCMV)
cd(savedir)
Run_sourceanalysis

%% Conn and network analyses (voxel and ROI levels)
if flag.connanalysis ==1
    if ~exist(fullfile(savedir,['wPLI', num2str(run), '_cond',num2str(condval),'.mat']), 'file') || ...
            ~exist(fullfile(savedir,['vs', num2str(run), '_cond',num2str(condval),'.mat']), 'file')
        Run_extractvirtualsensor %- Extract virtual sensors
        Run_connanalysis %- Conn analysis (wPLI)
        save(fullfile(savedir,['wPLI', num2str(run), '_cond',num2str(condval)]), 'source_conn1', 'par','vs_roi1','-v7.3');
        save(fullfile(savedir,['vs', num2str(run), '_cond',num2str(condval)]), 'vs','-v7.3');
    else
        load(fullfile(savedir,['wPLI', num2str(run), '_cond',num2str(condval),'.mat']));
        load(fullfile(savedir,['vs', num2str(run), '_cond',num2str(condval) ,'.mat']));
    end
end
Run_networkanalysis %- network analysis (network degrees)
Run_vis_netconn_analysis %- mapping

if flag.dicsanalysis == 1
    %% Source analysis, frequncy-domain (DICS-BF)
    Run_dicsanalysis
end
