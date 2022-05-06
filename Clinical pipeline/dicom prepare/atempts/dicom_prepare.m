




%% Read the DICOM files - With Nose - Low Res Dicom file 
%  helper-MRI (e.g. with nose)
% dicomfile = '/data/MEG/Clinical/MRI/xxx/DICOM/EXP00000/EXP0000';
cd('/MEG_data/Vahab/testdata/DICOM/Johnson_Ty_MEG_11-11-2019_DICOM/Johnson_Ty_MEG_11-11-2019_DICOM')
set_spm
dicom_left_dipole  = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');
dicom_right_dipole = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

% [pathstr, name] = fileparts(dicomfile);

set_ft
dipole_left = ft_read_mri(dicom_left_dipole);
dipole_right = ft_read_mri(dicom_right_dipole);
% ft_sourceplot([], mri_help);

%%
cd('/MEG_data/MRI_database/epilepsy/JOHNSON_Ty_Sag_SPGR/DICOM')

set_spm
dicomfile_MR = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

set_ft
mri_org = ft_read_mri(dicomfile_MR);
ft_sourceplot([], mri_org);

%%  Saving data as nii
savedir = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare';
cd(savedir)

nii_name = 'MRI_org.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_org);

nii_name = 'dipole_left.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, dipole_left);


nii_name = 'dipole_right.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, dipole_right);

%% Matlab function - test
% I = dicomread(dicomfile_MR);
% imtool(I,'DisplayRange',[])
% 
% info = dicominfo(dicomfile_MR);
% 
% cd(savedir)
% dic_name = 'test_dicom';
% dicomwrite(I,dic_name);
% 
% 
% I1 = dicomread(dic_name);
% info1 = dicominfo(dic_name);
% 
% I3 = int16(squeeze(mri_org.anatomy(:,:,1)));
% dic_name = 'test_dicom_1';
% dicomwrite(I3,dic_name);
% 
% imtool(I3,'DisplayRange',[])
% 
% info2 = dicominfo(dic_name);

%%
cd(savedir)

set_ft
dipole_fromBS = ft_read_mri('dipole_left_fromBS.nii');
T1 = ft_read_mri('T1.nii');
T1_mgz = ft_read_mri('/MEG_data/MRI_database/epilepsy/JOHNSON_Ty_Sag_SPGR/johnson_ty_sag_spgr_10_02_2019_FSRecon/mri/T1.mgz');
T1_org = ft_read_mri('/MEG_data/MRI_database/epilepsy/JOHNSON_Ty_Sag_SPGR/johnson_ty_sag_spgr_10_02_2019_FSRecon/mri/orig/001.mgz');

nii_name = 'T1_org.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, T1_org);

T1_org = ft_read_mri('T1_org.nii');
dipole_left_fromBS_dc = ft_read_mri('MRI_org.dcm');

dipole_test = ft_read_mri('/MEG_data/epilepsy/johnson_ty/191111/report/SpikesDicom/run2_1002_RightSpikes_dipoles_90gof_AntSensors.nii');

%%
set_ft
cd('/MEG_data/epilepsy/johnson_ty/191111/report/SpikesDicom')
T1 = ft_read_mri('T1.nii');
L_dip = ft_read_mri('run2_1001_LeftSpikes_dipoles_90gof_AntSensors.nii');

% cd('/MEG_data/epilepsy/johnson_ty/191111/report/SpikesDicom/DipolesDicom')
% T1 = ft_read_mri('T1_flirt.nii')

cd(savedir)
nii_name = 'L_dip.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, L_dip);

%% Reslice & save the transformation matrix to the anatomy_dir
% cfg                 = [];
% cfg.resolution      = 1;
% cfg.dim             = [256 256 150];
% dipole_fromBS_rs        = ft_volumereslice(cfg, dipole_fromBS);
% ft_sourceplot([], dipole_fromBS_rs);

%% SPM coreg, estimate and reslice
cd(savedir)

set_spm
nii_filename2 = 'T1_org.nii'; %256x256x150
nii_filename1 = 'T1.nii';      %256x256x256

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

% set_spm
% nii_filename2 = 'T1_org.nii'; %256x256x150
% nii_filename1 = 'L_dip.nii';  %256x256x256
% matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[nii_filename2,',1']};
% matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[nii_filename1,',1']};
% matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
% matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
% matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
% matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
% matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
% matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
% matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
% matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
% matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
% spm_jobman('run',matlabbatch);

spm_check_registration

%%
set_ft
rT1 = ft_read_mri('rT1.nii');
ft_sourceplot([], rT1);


cfg                 = [];
cfg.resolution      = 1;
cfg.dim             = [256 256 150];
sL_dip        = ft_volumereslice(cfg, L_dip);

rL_dip = ft_transform_geometry(T1_org.transform, sL_dip);
% rL_dip = ft_read_mri('L_dip.nii');
ft_sourceplot([], rL_dip);

%% resplice + oreg
% set_spm
% spm_coreg
% 
% spm_check_registration

%%
% metadata = dicominfo('CT-MONO2-16-ankle.dcm');
addpath('/MEG_data/Vahab/Github/MCW-MEGlab/FT/functions/External');
d = rdir('/MEG_data/MRI_database/epilepsy/JOHNSON_Ty_Sag_SPGR/DICOM/EXP00000/EXP*');

metadata_all = [];
for i=1:length(d)
    metadata = dicominfo(d(i).name);
    metadata_all{i}=metadata;
end

%%
cd(savedir)
cd ./dicom_test
for i=1:size(rT1.anatomy,3)
    I = int16(squeeze(rT1.anatomy(:,:,i)));
    if i < 11
        dicome_name = ['EXP000',num2str(i-1)];
    elseif i < 101
        dicome_name = ['EXP00',num2str(i-1)];
    else
        dicome_name = ['EXP0',num2str(i-1)];
    end
    %     dicomwrite(I,dicome_name, metadata_all{i});
    dicomwrite(I,dicome_name,metadata_all{i});
end
cd ..


cd(savedir)
cd ./dicom_leftdip
for i=1:size(sL_dip.anatomy,3)
    I = int16(squeeze(sL_dip.anatomy(:,:,i)));
    if i < 11
        dicome_name = ['EXP000',num2str(i-1)];
    elseif i < 101
        dicome_name = ['EXP00',num2str(i-1)];
    else
        dicome_name = ['EXP0',num2str(i-1)];
    end
    %     dicomwrite(I,dicome_name, metadata_all{i});
    dicomwrite(I,dicome_name,metadata_all{i});
end
cd ..

%%
set_spm
dicomfile_MR_test = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

set_ft
mri_test = ft_read_mri(dicomfile_MR_test);
ft_sourceplot([], mri_test);

nii_name = 'dicom_out_T1.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_test);


close all
set_spm
dicomfile_MR_test = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

set_ft
mri_test = ft_read_mri(dicomfile_MR_test);
ft_sourceplot([], mri_test);

nii_name = 'dicom_out_leftdip.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, mri_test);


