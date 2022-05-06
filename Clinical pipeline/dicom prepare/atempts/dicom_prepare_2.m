clear, clear, close all,

%%  Saving data as nii
savedir = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare';
cd(savedir)

%%
% cd('/MEG_data/epilepsy/johnson_ty/191111/report/SpikesDicom');
cd('/MEG_data/epilepsy/johnson_ty/191111/report/Spikes');

set_ft
T1_1 = ft_read_mri('T1.nii');
T1 = ft_read_mri('T1_BS.nii');
% T1_mgz = ft_read_mri('/MEG_data/MRI_database/epilepsy/JOHNSON_Ty_Sag_SPGR/johnson_ty_sag_spgr_10_02_2019_FSRecon/mri/T1.mgz');
L_dip = ft_read_mri('run2_1001_LeftSpikes_dipoles_90gof_AntSensors.nii');

%%
cd(savedir)
nii_name = 'L_dip.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, L_dip);
ft_sourceplot([], L_dip);

nii_name = 'T1.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, T1);
ft_sourceplot([], T1);

% nii_name = 'T1_mgz.nii';
% cfg                 = [];
% cfg.filename        = nii_name;
% cfg.filetype        = 'nifti';
% cfg.parameter       = 'anatomy';
% ft_volumewrite(cfg, T1_mgz);
% ft_sourceplot([], T1_mgz);

%% MR_org
cd('/MEG_data/MRI_database/epilepsy/JOHNSON_Ty_Sag_SPGR/DICOM')

set_spm
dicomfile_MR = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

set_ft
dicom_MR = ft_read_mri(dicomfile_MR);
ft_sourceplot([], dicom_MR);

%%
cd(savedir)
nii_name = 'dicom_MR.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, dicom_MR);

%%
% T1_flirt = ft_read_mri('T1_flirt.nii');

%% SPM coreg, estimate and reslice
cd(savedir)

set_spm
nii_filename2 = 'dicom_MR.nii'; % 256x256x150
nii_filename1 = 'T1.nii';       % 256x256x256

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
% nii_filename2 = 'dicom_MR.nii'; % 256x256x150
% nii_filename1 = 'L_dip.nii';       % 256x256x256
% 
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

% spm_check_registration

%% Reslise
nii_filename2 = 'dicom_MR.nii'; % 256x256x150
nii_filename1 = 'L_dip.nii';       % 256x256x256

flag = [];
flag.mask = 0;
flag.mean = 0;
flag.interp = 4;
flag.which = 1;
flag.wrap = [0 0 0];
flag.prefix = 'r';
P = {nii_filename1;nii_filename2};
spm_reslice(P, flag)

%%
set_ft
rT1 = ft_read_mri('rT1.nii'); 
rL_dip = ft_read_mri('rL_dip.nii');

ft_sourceplot([], rT1);
ft_sourceplot([], rL_dip);

%%
% cfg                 = [];
% cfg.resolution      = 1;
% cfg.dim             = [256 256 150];
% sL_dip        = ft_volumereslice(cfg, L_dip);

% rL_dip = ft_transform_geometry(rT1.transform/sL_dip.transform, sL_dip);
% rtL_dip = ft_transform_geometry(rT1.transform, rL_dip);
rtL_dip = ft_transform_geometry(rT1.transform/rL_dip.transform, rL_dip);

% rL_dip = ft_read_mri('L_dip.nii');
ft_sourceplot([], rL_dip);

%%
[a,b] = fileparts(dicomfile_MR);
addpath('/MEG_data/Vahab/Github/MCW-MEGlab/FT/functions/External');
d = rdir(fullfile(a,'EXP*'));

metadata_all = [];
for i=1:length(d)
    metadata = dicominfo(d(i).name);
    metadata_all{i}=metadata;
end

%%
outd = fullfile(savedir,'dicom','T1'); if exist(outd, 'file') == 0, mkdir(outd); end

cd(savedir)
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
    dicomwrite(I,fullfile(outd,dicome_name), metadata_all{i});
end


outd = fullfile(savedir,'dicom','Ldip'); if exist(outd, 'file') == 0, mkdir(outd); end
for i=1:size(rL_dip.anatomy,3)
    I = int16(squeeze(rL_dip.anatomy(:,:,i)));
    if i < 11
        dicome_name = ['EXP000',num2str(i-1)];
    elseif i < 101
        dicome_name = ['EXP00',num2str(i-1)];
    else
        dicome_name = ['EXP0',num2str(i-1)];
    end
    %     dicomwrite(I,dicome_name, metadata_all{i});
    dicomwrite(I,fullfile(outd,dicome_name),metadata_all{i});
end

%%
set_spm
dicomfile_T1_check = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

set_ft
T1_check = ft_read_mri(dicomfile_T1_check);
ft_sourceplot([], T1_check);

nii_name = 'dicom_out_T1.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, T1_check);

close all
set_spm
dicomfile_Ldip_check = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

set_ft
ldip_check = ft_read_mri(dicomfile_Ldip_check);
ft_sourceplot([], ldip_check);

nii_name = 'dicom_out_ldip.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, ldip_check);

%% Validation test




