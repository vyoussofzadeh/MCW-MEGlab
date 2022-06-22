%% The bgros_pilot

% MEG (pre)-processing pipeline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 05/24/2022

clear; clc, close('all'); warning off

%% Flag
flag = [];
flag.plottriggers = 1;
flag.anatomy_check = 1;
flag.inspect = 0;

%% Paths
restoredefaultpath
script_path = '/group/bgross/work/CIDMEG/analysis/FT pipelines';
addpath(genpath(script_path));

%- Input dir
indir = '/group/bgross/work/CIDMEG/RawData/cidmeg_1/220517/tsss';
%- Output dir
outdir = '/group/bgross/work/CIDMEG/analysis/process';

%
ft_path = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_20190419');
addpath(ft_path);
ft_defaults

allpath = [];
allpath.ft_path18 = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_041718'); % needed for IC plotting
allpath.ft_path = ft_path;
allpath.script_path = script_path;

%- atlas
atlas_path = fullfile(ft_path,'template','atlas');
atlas = ft_read_atlas(fullfile(atlas_path,'aal/ROI_MNI_V4.nii'));

template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %

%% anatomy settings
mridir = '/group/bgross/work/CIDMEG/analysis/anatomy/Sub001/T1';
mrifile = 'sub-CIDFMRI001_ses-1_acq-mprage_T1w.nii.gz';

%% Reading data, Layout & sensor location
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);

%% Enter subject ID, run, and condition info
subj = 'pilot1';
run = input('enter run of intersts, e.g, 1: ');
switch run
    case {6, 7}
        run = '6_7';
end
condval = input('enter condition of intersts, e.g, 1, 2, 4, 8: ');
datafilename = ['Task_run', num2str(run), '_raw.fif'];
datafile = fullfile(indir, datafilename);

subdir = fullfile(outdir, subj,'preprocess'); % output dir
if exist(subdir, 'file') == 0, mkdir(subdir); end
cd(subdir)

if ~exist(fullfile(subdir,['data_run', num2str(run), '_cond',num2str(condval), '.mat']), 'file')
    
    %%
    sens = ft_read_sens(datafile);
    sens = ft_convert_units(sens,'mm');
    
    %% Reading trigger info
    cfg = [];
    cfg.plot = 1;
    [detResp,detTrig] = do_readtriggers(cfg, datafile);
    
    %- THIS IS OPTIOPNAL
    evt = ft_read_event(datafile);
    
    %% Preprocessing, band-pass filtering
    cfg = [];
    cfg.datafile = datafile;
    cfg.prestim = 1;
    cfg.poststim = 2;
    % cfg.eventvalue = [1, 2, 4, 8]; % this can be changed based on the event type (see event var)
    cfg.eventvalue = condval; % this can be changed based on the event type (see event var)
    cfg.hpfreq  = 0.1;
    cfg.lpfreq = 100;
    cfg.dftfreq = 60; % or [60 120 180];
    f_data = do_ft_preprocess(cfg);
    
    %% Preprocessing, notch filtering
    cfg = [];
    cfg.subj = subj;
    fn_data = do_notch(cfg, f_data);
    
    %% Preprocessing, ICA
    cfg = [];
    cfg.lay = lay;
    cfg.subj = num2str(run);
    cfg.n = 20;
    cfg.allpath = allpath;
    cfg.savefig = 1;
    ica_data = do_ica(cfg, fn_data);
    close all
    
    %% Preprocessing, rejecting bad trials
    cfg = [];
    cfg.pflag = 1; % yes:1, No:2
    cfg.saveflag = 0; % yes:1, No:2
    cfg.savepath = [];
    cfg.latency = [ica_data.time{1}(1),ica_data.time{1}(end)];
    cfg.rejectpercentage = .95;
    cfg.method = 'auto'; % 'manual'
    [data_clean,report] = do_rejectdata(cfg, ica_data);
    
    %% Save preprocessed data
    disp('saving data ...')
    trl = []; trl.all = f_data.trial; trl.cleaned = data_clean.trial; trl.report.btrl = report.btrl;
    cd(outdir)
    save(fullfile(subdir,['data_run', num2str(run), '_cond',num2str(condval)]), 'data_clean', 'trl', '-v7.3');
    
else
    cd(outdir)
    disp('data already preprocessed!')
%     load(fullfile(subdir,['data_run', num2str(run), '_cond',num2str(condval), '.mat']))
end

close all,
%% Inspecting data
if flag.inspect ==1
    
    % Freq analysis
    cfg = [];
    cfg.savefile = [];
    cfg.saveflag = 0;
    cfg.foilim = [2 50];
    cfg.plotflag  = 1;
    cfg.tapsmofrq       = 4;
    cfg.taper    = 'hanning';
    do_fft(cfg, data_clean);
    
    % Time-Freq analysis
    cfg = [];
    cfg.layout = lay;
    cfg.subj = subj;
    cfg.baseline = [-0.3,0];
    cfg.fmax = 50;
    cfg.title = 'Time-Freq';
    [t_max,f_max, tfr]  = do_tfr_analysis(cfg, data_clean);
    
    % Grand Mean
    a_data = do_ave(data_clean);
    cfg = [];
    cfg.savefile = [];
    cfg.saveflag = 2;
    cfg.lay  = lay;
    do_ave_plot(cfg, a_data);
    
end