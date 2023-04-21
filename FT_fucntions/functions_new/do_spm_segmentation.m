function mri_brain = do_spm_segmentation(cfg_main)

nii_file = cfg_main.nii_file;
spm_path =  cfg_main.spm_path;

%%
% Add SPM12 to MATLAB path
% spm_path = 'path/to/spm12'; % Replace with your SPM12 path
addpath(spm_path);
spm('defaults', 'PET');

% Load NIfTI file
% nii_file = 'path/to/your/file.nii'; % Replace with your NIfTI file path
vol = spm_vol(nii_file);
vol_data = spm_read_vols(vol);

matlabbatch = [];
matlabbatch{1}.spm.spatial.preproc.channel.vols = {nii_file}; % Input NIfTI file
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,1')}; % TPM tissue probability map for gray matter
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,2')}; % TPM tissue probability map for gray matter
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,3')}; % TPM tissue probability map for gray matter
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,4')}; % TPM tissue probability map for gray matter
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,5')}; % TPM tissue probability map for gray matter
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm_path, 'tpm', 'TPM.nii,6')}; % TPM tissue probability map for gray matter
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];
% Run the segmentation job
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

% Check for successful completion
fprintf('Segmentation completed successfully.\n');

% Remove SPM12 from MATLAB path
rmpath(spm_path);

%%
mri_gray = ft_read_mri(['c1', nii_file]); ft_sourceplot([], mri_gray);
mri_csf = ft_read_mri(['c2', nii_file]); ft_sourceplot([], mri_csf);
mri_white = ft_read_mri(['c3', nii_file]); ft_sourceplot([], mri_white);

%%
segmented = [];
segmented.gray = mri_gray.anatomy;
segmented.csf = mri_csf.anatomy;
segmented.white = mri_white.anatomy;

% create the brain from the tpm
fprintf('creating brainmask ... using the summation of gray, white and csf tpms\n');
brain = segmented.gray + segmented.white + segmented.csf;
brainmask = brain>0;
segmented.brain = brainmask;

mri_brain = mri_gray; mri_brain.brain = segmented.brain;
% 
% 
% cfg = [];
% cfg.method = 'singleshell';
% cfg.spmversion = 'spm12';
% individual_headmodel = ft_prepare_headmodel(cfg, mri_brain1);
% 
% figure;
% ft_plot_mesh(individual_headmodel.bnd, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;