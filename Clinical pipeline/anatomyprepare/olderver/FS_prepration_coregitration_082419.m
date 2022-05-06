clear, clear, close all,

%%
cd_org = '/MEG_data/Vahab/Shared Scripts/Freesurfer_anat_prepration/shared';
addpath(cd_org)
path_tools = '/usr/local/MATLAB_Tools';

%%
nii_filename1 = '/MEG_data/MRI_database/epilepsy/HAINES_Robert_Nose/mni_resliced.mgz';
nii_filename2 = '/MEG_data/MRI_database/epilepsy/HAINES_Robert_1/mni_resliced.mgz';

%%
set_ft
mri1 = ft_read_mri(nii_filename1);
ft_sourceplot([], mri1);
mri2 = ft_read_mri(nii_filename2);
ft_sourceplot([], mri2);

%% Save the resliced mni-transformed mri image
nii_filename1 = 'mni_resliced1.nii';
cfg                 = [];
cfg.filename        = nii_filename1;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri1);

nii_filename2 = 'mni_resliced2.nii';
cfg                 = [];
cfg.filename        = nii_filename2;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri2);

%% Co-reg, estimate and reslice
% restoredefaultpath
% addpath(cd_org)
% spm_path = fullfile(path_tools,'/spm12');
% addpath(genpath(spm_path))
% spm_get_defaults
set_spm

%% SPM coreg, estimate and reslice
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[nii_filename2,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[nii_filename1,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);
spm_coreg

spm_check_registration
% spm_image

%% ave of 
spm_mean_ui

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

%% bias correction on high-res MRI, no nose

% set_spm8
biasfield = spm_bias_estimate('mni_resliced2.nii');
spm_bias_apply('mni_resliced2.nii', biasfield);

