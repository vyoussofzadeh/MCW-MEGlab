
clear, clear, clc, close all,

% Author, Vahab Youssof Zadeh, 2021
% update: 09/27/22,  scp copy files to Squiggles was added

%%
cd_org = '/MEG_data/MCW_pipeline/Anatomyprepare';
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
cd(pathstr)

%%
set_ft
mri = ft_read_mri(dicomfile);

%%
ft_sourceplot([], mri);
mri = ft_convert_units(mri, 'mm');

%% Filename for saving
clc
cd(pathstr)
cd ..
cd ..
disp(['1 = suggesting:', pwd])
disp( '2 = other')
indir_ask = input('Enter saveing dir:');

switch indir_ask
    case 1
        indir = pwd;
    case 2
        indir = input('Enter saveing dir:','s');
end
cd(indir)
nii_filename = 'mri.nii';

nii_savepath = 'nii';
if exist(nii_savepath, 'file') == 0, mkdir(nii_savepath), end

%% Save the resliced mni-transformed mri image
cfg                 = [];
cfg.filename        = fullfile(nii_savepath, nii_filename);
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri);

%%
cd(nii_savepath)
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
cd ..
cfg                 = [];
cfg.filename        = fullfile(nii_savepath,'mri_resliced.nii');
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_resliced);

%% Save the transformation matrix
% transform_vox2mni   = mri_resliced.transform;
% filename_vox2mni    = 'transform_vox2mni';
% save(filename_vox2mni, 'transform_vox2mni');

%% Bias correction
disp('========');
disp('1: Yes');
disp('2: No');
bsask = input('Bias correction?');
if bsask ==1
    set_spm
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

%% Bias correction
cd(indir)
clc
disp('========');
disp('1: Yes');
disp('2: No');
fsask = input('Run FS?');
% subn = input('enter subject name:');

[pathstr, name, ext] = fileparts (pwd);

disp(['1 = suggesting:', name])
disp( '2 = other')
ask_subn = input('Enter saveing dir:');

switch ask_subn
    case 1
        subn = name;
    case 2
        subn = input('enter subject name:');
end

% savepath = 'FS';
% if exist(savepath, 'file') == 0, mkdir(savepath), end
% FS_savedir = fullfile(pathstr, subn, savepath);

if fsask ==1
    clc, close all,
    disp('run this in command ...'),
    disp(['cd ', fullfile(indir, nii_savepath)]);
%     [a, name] = fileparts(pwd);
    disp(['recon-all -s ', subn, ' -i mri_resliced.nii -all'])
%     disp(['recon-all -s ', FS_savedir, ' -i mri_resliced.nii -all'])
%     command = ['recon-all -s ', name, ' -i mri_resliced.nii -all'];
%     system(command)
end


%% Check FS
disp('1: Yes');
disp('2: No');
fschkask = input('Check FS?');
if fschkask ==1
    % /MEG_data/MRI_database/epilepsy/RAPEY_Ward_Jacob/FSrecon_110921
    % cd /MEG_data/MRI_database/epilepsy/RAPEY_Ward_Jacob/
    % tkmedit FSrecon_110921 T1.mgz -surfs
    
end

%% Copy FS
disp('==========')
disp('1: Yes');
disp('2: No');
catask = input('Copy nii to Sqioggles for Cat analysis?');
if catask
    % - copy data between servers
    command = (['scp -r ', fullfile(pwd,'nii'), ' vyoussofzadeh@squiggles.rcc.mcw.edu:/data/MEG/Vahab/Data_clinical/CAT_analysis/', name]);
    system(command)
    disp('copied to /data/MEG/Vahab/Data_clinical/CAT_analysis/')
end


