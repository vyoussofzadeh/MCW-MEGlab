
clear, clear, close all,

% ===========================================
% Freesurfer segmenation prepration
% Author, Vahab Youssof Zadeh,
% Last update: 07/29/2020
% ===========================================

%%
disp('1: For CT-MR coreg, grid/SEEG processing')
disp('2: For MRI segmentation (MEG processing)')
in_sel = input(':');

switch in_sel
    case 1
        indir = '/MEG_data/JEFF/GRID_Processing/';
    case 2
        indir = '/MEG_data/MRI_database/epilepsy/';
end

cd(indir)
[subjdir] = uigetdir;
cd(subjdir)

%%
cd_org = '/MEG_data/Vahab/Shared Scripts/Freesurfer_anat_prepration/shared';
addpath(cd_org)
path_tools = '/usr/local/MATLAB_Tools';

%% org-MRI (e.g. High Res, no-nose)
nii_filename_org = fullfile(subjdir,'mri/T1.mgz'); 

set_ft
mri_org = ft_read_mri(nii_filename_org);
ft_sourceplot([], mri_org);

%% Read the DICOM files - With Nose - Low Res Dicom file 
%  helper-MRI (e.g. with nose)
% dicomfile = '/data/MEG/Clinical/MRI/xxx/DICOM/EXP00000/EXP0000';

set_spm
dicomfile = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');
[pathstr, name] = fileparts(dicomfile);

set_ft
mri_help = ft_read_mri(dicomfile);
ft_sourceplot([], mri_help);

%%
cd ..
savedir = 'MRCoreg_spm';
if exist(savedir, 'file') == 0, mkdir(savedir); end
cd(savedir)

%% Save the resliced mni-transformed mri image
nii_mri_org = 'MRI_org.nii';
cfg                 = [];
cfg.filename        = nii_mri_org;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_org);

nii_mri_helper = 'MRI_helper.nii';
cfg                 = [];
cfg.filename        = nii_mri_helper;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_help);

%% spm CoReg, estimate and reslice
set_spm
spm_coreg

%%
spm_check_registration

%% Ave of MRI_org and MRI_helper
spm_mean_ui

%%
spm_check_registration

%%
set_ft
addpath('/usr/local/MATLAB_Tools/fieldtrip_20190419/external/spm8')

%% spm8 bias correction
% set_spm8
biasfield = spm_bias_estimate('mean.nii');
spm_bias_apply('mean.nii', biasfield);

%%
mean_nii = ft_read_mri('mmean.nii');
nii_filename_final = 'T1.mgz';
cfg                 = [];
cfg.filename        = nii_filename_final;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mean_nii);
