     
close all
clc

%%
addpath('/usr/local/MATLAB_Tools/mne')
addpath('/MEG_data/MCW_pipeline/Preprocess/func')

%%
cd('/MEG_data/epilepsy/belcourt_ashley/brainstorm_db/data')
rawFile  = '/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run02_spont_supine/Run02_spont_supine_raw_t_sss_ecgClean_raw_BAK.fif'; 
outFile  = '/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run02_spont_supine/Run02_spont_supine_raw_t_sss_ecgClean_raw_fixed.fif';
Ch = load('./belcourt_ashley/@rawRun02_spont_supine_raw_t_sss_ecgClean_raw_BAK/channel_vectorview306_acc1.mat');


rawFile  = '/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run03_spont_supine/Run03_spont_supine_raw_t_sss_ecgClean_raw_BAK.fif'; 
outFile  = '/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run03_spont_supine/Run03_spont_supine_raw_t_sss_ecgClean_raw_fixed.fif';
Ch = load('./belcourt_ashley/@rawRun03_spont_supine_raw_t_sss_ecgClean_raw_BAK/channel_vectorview306_acc1.mat');

rawFile  = '/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run04_spont_supine/Run04_spont_supine_raw_t_sss_ecgClean_raw_BAK.fif'; 
outFile  = '/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run04_spont_supine/Run04_spont_supine_raw_t_sss_ecgClean_raw_fixed.fif';
Ch = load('./belcourt_ashley/@rawRun04_spont_supine_raw_t_sss_ecgClean_raw_BAK/channel_vectorview306_acc1.mat');

rawFile  = '/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run05_spont_supine/Run05_spont_supine_raw_t_sss_ecgClean_raw_BAK.fif'; 
outFile  = '/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run05_spont_supine/Run05_spont_supine_raw_t_sss_ecgClean_raw_fixed.fif';
Ch = load('./belcourt_ashley/@rawRun05_spont_supine_raw_t_sss_ecgClean_raw_BAK/channel_vectorview306_acc1.mat');

rawFile  = '/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run06_spont_supine/Run06_spont_supine_raw_t_sss_ecgClean_raw_BAK.fif'; 
outFile  = '/MEG_data/epilepsy/belcourt_ashley/250612/sss/Run06_spont_supine/Run06_spont_supine_raw_t_sss_ecgClean_raw_fixed.fif';
Ch = load('./belcourt_ashley/@rawRun06_spont_supine_raw_t_sss_ecgClean_raw_BAK/channel_vectorview306_acc1.mat');

%%
% 1) load the Brainstorm channel file you showed
T_dev2head_mm = Ch.TransfMeg{2};          % 4×4, millimetres

% --------------------------------------------------------------
% 2. build +90° rotation around +Z  (clockwise when viewed from above)
%    Rz(+90) converts Brainstorm frame into Neuromag head frame
% --------------------------------------------------------------
Rz = [ 0 -1  0  0;     %   [  0  -1   0   0 ]
      +1  0  0  0;     %   [ +1   0   0   0 ]
       0  0  1  0;     %   [  0   0   1   0 ]
       0  0  0  1];    %   [  0   0   0   1 ]

% --------------------------------------------------------------
% 3. apply rotation  ➜ new Device→Head (still mm)
% --------------------------------------------------------------
T_fix_mm = Rz * T_dev2head_mm;

% --------------------------------------------------------------
% 4. convert translations to **metres** for Neuromag FIF
% --------------------------------------------------------------
T_fix_m          = T_fix_mm;
T_fix_m(1:3,4)   = T_fix_m(1:3,4) / 1000;   % mm → m

%%
restoredefaultpath
path_tools = '/usr/local/MATLAB_Tools';
ft_path = fullfile(path_tools,'/fieldtrip_2022');

addpath(ft_path);
ft_defaults

hdr   = ft_read_header(rawFile,'headerformat','neuromag_mne');
hsp   = ft_convert_units(ft_read_headshape(rawFile),'cm');   % grey points
grad  = ft_convert_units(ft_transform_sens(T_fix_m,hdr.grad),'cm'); % blue

figure;
ft_plot_headshape(hsp,'vertexcolor',[.7 .7 .7],'vertexsize',4); hold on
ft_plot_sens(grad,'style','*b');
view([90 0]); axis equal; title('After +90° Z-rotation');

%%
addpath('/usr/local/MATLAB_Tools/mne')
addpath('/MEG_data/MCW_pipeline/Preprocess/func')

cfg = struct('infile',   rawFile, ...
             'outfile',  outFile, ...
             'T',        T_fix_m, ...     % final 4×4 (metres)
             'blocksec', 10);

do_mne_write_with_transform(cfg);

% quick confirmation
info = fiff_read_meas_info(outFile);
disp(info.dev_head_t.trans);              % should equal T_fix_m



