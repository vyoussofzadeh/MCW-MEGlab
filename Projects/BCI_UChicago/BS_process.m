%% BCI data, Medical College of Wisconsin - Uni Chicago

% Script: BS Process (preprocessing, source analysis)
% Project: BCI
% Writtern by: Vahab YoussofZadeh
% Update: 03/26/2023

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
cd('/MEG_data/Research_studies/BCI_Uni_of_Chicago')

indir = '/MEG_data/Research_studies/BCI_Uni_of_Chicago/bci03_bci03';
outdir = '/MEG_data/Research_studies/BCI_Uni_of_Chicago/Process';

ft_path = '/usr/local/MATLAB_Tools/fieldtrip_20190419';
addpath(ft_path); ft_defaults

%%
ft_func = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/FT_fucntions';
addpath(fullfile(ft_func, '/External/brewermap'));
addpath(fullfile(ft_func, '/functions_new'));

%%