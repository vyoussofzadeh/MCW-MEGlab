clear, close all,
%% 
restoredefaultpath
addpath(pwd)
%- Adding path
path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
% cfg_init.path_tools = '/home/vyoussof/Desktop/Vahab/Scripts/tools';
% [allpath, ~] = vy_init(cfg_init);

%- Adding path
cfg_init = [];
% cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
% cfg_init.path_tools = '/home/vyoussof/Desktop/Vahab/Scripts/tools';
% cfg_init.path_tools = '/usr/local/MATLAB_Tools';
% ft_path = fullfile(path_tools,'/fieldtrip');
ft_path = fullfile(path_tools,'/ft_packages/fieldtrip_master');

addpath(ft_path);
ft_defaults

% spm_path = fullfile(cfg.path_tools,'/spm12');
spm_path = fullfile(path_tools,'SPM/spm12');
addpath(spm_path);

%% Read the DICOM files
% dicomfile = '/data/MEG/Clinical/MRI/xxx/DICOM/EXP00000/EXP0000';
% [f, p] = uigetfile('*'); dicomfile = fullfile(p, f);
dicomfile = spm_select(1,'.*','Select one dicome file');
[pathstr, name] = fileparts(dicomfile);
%%
mri = ft_read_mri(dicomfile);
%%
ft_sourceplot([], mri);
mri = ft_convert_units(mri, 'mm');

%% Filename for saving
cd(pathstr)
cd ..
cd ..
nii_filename = 'mri.nii';

%% Save the resliced mni-transformed mri image
cfg                 = [];
cfg.filename        = nii_filename;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri);

%%
% addpath(genpath(allpath.spm_path))
spm_get_defaults
setorigin_center(nii_filename);

%%
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

%% Filename for saving
mgz_filename = 'mni_resliced.mgz';

%% Save the resliced mni-transformed mri image
cfg                 = [];
cfg.filename        = mgz_filename;
cfg.filetype        = 'mgz';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_resliced);

%% Save the transformation matrix
transform_vox2mni   = mri_resliced.transform;
filename_vox2mni    = 'transform_vox2mni';
save(filename_vox2mni, 'transform_vox2mni');

%%
