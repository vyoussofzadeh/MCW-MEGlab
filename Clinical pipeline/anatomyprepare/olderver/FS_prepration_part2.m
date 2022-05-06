clear, clear, close all,

%%
% cd_org = '/MEG_data/Vahab/Shared Scripts/Freesurfer_anat_prepration/shared';
cd_org = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/anatomyprepare';
addpath(cd_org)
path_tools = '/usr/local/MATLAB_Tools';

%%
set_spm
niifile = spm_select(1,'.nii','Select nii');

%%
set_ft

%% Bias correction
% set_spm8
biasfield = spm_bias_estimate(niifile);
spm_bias_apply('mni_resliced.nii', biasfield);

%%
mni_resliced = ft_read_mri('mni_resliced.nii');
ft_sourceplot([], mni_resliced); title('before BS correction')
mri_biascorrected = ft_read_mri('mmni_resliced.nii');
ft_sourceplot([], mri_biascorrected); title('after BS correction')

