%% BCI dataset, Ulster University & Medical College of Wisconsin

% Script: setup paths
% Project: BCI_nonstationarity
% Writtern by: Vahab Youssof Zadeh
% Update: 11/03/2022

%%
restoredefaultpath
script_path = '/data/MEG/Vahab/Github/MCW_MEGlab/Projects/BCI/non_stationary';
addpath(genpath(script_path));

%- Input dir
indir = '/data/MEG/Vahab/Github/MCW_MEGlab/Projects/BCI/Data/MEG_mat';
%- Output dir
outdir = '/data/MEG/Vahab/Github/MCW_MEGlab/Projects/BCI/Process/ft_process';

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
