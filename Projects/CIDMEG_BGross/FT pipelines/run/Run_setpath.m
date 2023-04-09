restoredefaultpath
script_path = '/group/bgross/work/CIDMEG/analysis/FT pipelines';
addpath(genpath(script_path));

%- Input dir
indir = '/group/bgross/work/CIDMEG/RawData/cidmeg_1/220517/tsss';
%- Output dir
outdir = '/group/bgross/work/CIDMEG/analysis/process';

%
% ft_path = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_20190419');
ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
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