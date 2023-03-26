%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: Process (extract virtual sensrors at voxel and roi AAL atlas levels)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 11/16/2022

%%
%- Input dir
indir = '/group/jbinder/ECP/MEG/MEG_Work';
%- Output dir
outdir = '/data/MEG/Research/ECP/Semantic_Decision/process';

%
% ft_path = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_20190419');
ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';

allpath = [];
allpath.ft_path18 = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_041718'); % needed for IC plotting
allpath.ft_path = ft_path;

%%
cd('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/brewermap')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/Colormaps-from-MatPlotLib2.0')
addpath(ft_path); ft_defaults
