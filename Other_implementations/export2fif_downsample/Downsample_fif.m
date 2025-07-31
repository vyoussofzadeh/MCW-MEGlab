close all
clc

%%
% addpath('/usr/local/MATLAB_Tools/mne')
% addpath('/MEG_data/MCW_pipeline/Preprocess/func')
% 
% 
% cd('/MEG_data/epilepsy/belcourt_ashley/brainstorm_db/data')
% rawFile  = '/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run02_spont_supine/Run02_spont_supine_raw_t_sss_ecgClean_raw_BAK.fif'; 

ft_path = '/opt/matlab_toolboxes/ft_packages/fieldtrip_latest';
ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
addpath(ft_path);
ft_defaults


cd('/data/MEG/Research/SpikeDectection/tsss_ecgClean_Data')


cfg          = [];
cfg.dataset  = 'abt_marisa_Run02_spont_supine_raw_t_sss_ecgClean_raw.fif';
cfg.demean   = 'yes';      % optional
data         = ft_preprocessing(cfg);

cfg          = [];
cfg.resamplefs = 200;
data_ds      = ft_resampledata(cfg, data);

% Now hand data_ds to Brainstorm or MNE-MATLAB to write FIF

cfg = [];
cfg.infil = 'abt_marisa_Run02_spont_supine_raw_t_sss_ecgClean_raw.fif';
cfg.outfile = 'abt_marisa_Run02_spont_supine_raw_t_sss_ecgClean_raw_DS.fif';
cfg.cln_data = data_ds;
do_mne_ex_read_write_raw(cfg)

%%
cd('/data/MEG/Research/SpikeDectection/tsss_ecgClean_Data')
addpath('/opt/mne')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Other_implementations/export2fif_downsample')

cfg = [];
cfg.infile = 'abt_marisa_Run02_spont_supine_raw_t_sss_ecgClean_raw.fif';
cfg.outfile = 'abt_marisa_Run02_spont_supine_raw_t_sss_ecgClean_raw_DS.fif';
cfg.sfreq_new = 200;
cfg.blocksec   = inf;           % optional; default is 10 s
do_mne_write_downsample_clean(cfg)

%%
filename = cfg.outfile;
hdr = ft_read_header(filename);

filename = cfg.infile;
hdr = ft_read_header(filename);


[info, ~] = fiff_read_meas_info(filename);
fprintf('Fs = %.1f Hz\n', info.sfreq);
