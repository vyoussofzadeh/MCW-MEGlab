clear, clear, clc, close all,

%%
% cd_org = '/MEG_data/Vahab/Shared Scripts/Freesurfer_anat_prepration/shared';
cd_org = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/anatomyprepare';
addpath(cd_org)
path_tools = '/usr/local/MATLAB_Tools';

%% Read the DICOM files
% dicomfile = '/data/MEG/Clinical/MRI/xxx/DICOM/EXP00000/EXP0000';
% [f, p] = uigetfile('*'); dicomfile = fullfile(p, f);
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
set_spm
dicomfile = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');
[pathstr, name] = fileparts(dicomfile);

%%
set_ft
mri = ft_read_mri(dicomfile);

%%
ft_sourceplot([], mri);
mri = ft_convert_units(mri, 'mm');

%% Filename for saving
cd(pathstr)
indir = input('Enter saveing dir:');
cd(indir)
nii_filename = 'mri.nii';

%% Save the resliced mni-transformed mri image
cfg                 = [];
cfg.filename        = nii_filename;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri);

%%
set_spm
setorigin_center(nii_filename);

%%
set_ft
mri1 = ft_read_mri(nii_filename);
ft_sourceplot([], mri1);

%% Reslice & save the transformation matrix to the anatomy_dir
cfg                 = [];
cfg.resolution      = 1;
cfg.dim             = [256 256 256];
mri_resliced        = ft_volumereslice(cfg, mri1);
ft_sourceplot([], mri_resliced);

%%
% cfg          = [];
% cfg.method   = 'interactive';
% cfg.coordsys = 'spm';
% mri_mni     = ft_volumerealign(cfg, mri_resliced);

%% Filename for saving (.mgz)
% mgz_filename = 'mni_resliced.mgz';
%
% cfg                 = [];
% cfg.filename        = mgz_filename;
% cfg.filetype        = 'mgz';
% cfg.parameter       = 'anatomy';
% ft_volumewrite(cfg, mri_resliced);

%% Filename for saving (.nii)
nii_filename = 'mri_resliced.nii';
cfg                 = [];
cfg.filename        = nii_filename;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_resliced);

%% Save the transformation matrix
transform_vox2mni   = mri_resliced.transform;
filename_vox2mni    = 'transform_vox2mni';
save(filename_vox2mni, 'transform_vox2mni');

%% Bias correction
disp('1: Yes');
disp('2: No');
bsask = input('Bias correction?');

if bsask ==1
    
%     % - SPM-8
%     biasfield = spm_bias_estimate(nii_filename);
%     spm_bias_apply('mni_resliced.nii', biasfield);
%     
%     %
%     set_ft
%     mni_resliced = ft_read_mri('mni_resliced.nii');
%     ft_sourceplot([], mni_resliced); title('before BS correction')
%     mri_biascorrected = ft_read_mri('mmni_resliced.nii');
%     ft_sourceplot([], mri_biascorrected); title('after BS correction')
    
    %- SPM-12
    % https://layerfmri.com/2017/12/21/bias-field-correction/
    
    %-----------------------------------------------------------------------
    % Job saved on 09-Nov-2017 14:54:58 by cfg_util (rev $Rev: 6460 $)
    % spm SPM - SPM12 (6906)
    % cfg_basicio BasicIO - Unknown
    %-----------------------------------------------------------------------
    addpath('/usr/local/MATLAB_Tools/spm12')
    spm_defaults
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {'./mri_resliced.nii,1'};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 20;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
    spm_jobman('run',matlabbatch);
end
