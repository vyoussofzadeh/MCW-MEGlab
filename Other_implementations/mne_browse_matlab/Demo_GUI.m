restoredefaultpath
path_tools = '/usr/local/MATLAB_Tools';
ft_path = fullfile(path_tools,'/fieldtrip_2022');

addpath(ft_path);
ft_defaults

addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Other_implementations/mne_browse_matlab')

% cd '/MEG_data/epilepsy/magolan_thomas/250521/sss/Run02_spont_supine'
% mne_browser_gui('Run02_spont_supine_ic_raw_t_sss.fif')

cd('/MEG_data/epilepsy/aranda_ana/220224/sss/Run07_spont_supine')
mne_browser_gui('Run07_spont_supine_raw_t_sss_ecgClean_raw.fif')

cd('/MEG_data/epilepsy/aranda_ana/220224/sss/Run07_spont_supine')
mne_browser_gui_V1_53('Run07_spont_supine_raw_t_sss_ecgClean_raw.fif')


cd('/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run02_spont_supine')
mne_browser_gui('Run02_spont_supine_raw_sss_ecgClean_raw.fif')

cd('/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run02_spont_supine')
mne_browser_gui('Run02_spont_supine_raw_t_sss_ecgClean_raw.fif')


clc, clear, close all
cd('/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run02_spont_supine')
mne_browser_gui_V1_56('Run02_spont_supine_raw_t_sss_ecgClean_raw.fif')

system('mbrowse /MEG_data/epilepsy/belcourt_ashley/250612/sss/Run02_spont_supine/Run02_spont_supine_raw_t_sss_ecgClean_raw.fif')


clc, clear, close all
cd('/MEG_data/epilepsy/oates_latasha/250709/sss/Run02_spont_supine');
mne_browser_gui_V1_57('Run02_spont_supine_raw_t_sss_ecgClean_raw.fif')

system('mbrowse /MEG_data/epilepsy/oates_latasha/250709/sss/Run02_spont_supine/Run02_spont_supine_raw_t_sss_ecgClean_raw.fif')
