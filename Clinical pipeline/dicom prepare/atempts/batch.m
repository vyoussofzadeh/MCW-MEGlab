% List of open inputs
nrun = X; % enter the number of runs here
jobfile = {'/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare/batch_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});



matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {'/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare/dicom_MR.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {'/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare/T1.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {'/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare/L_dip.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
