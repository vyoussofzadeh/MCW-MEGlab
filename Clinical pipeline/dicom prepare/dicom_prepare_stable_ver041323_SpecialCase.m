
clc, clear, close all,
% set(0,'DefaultFigureWindowStyle','normal')

%%
cd_org = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Clinical pipeline/dicom prepare';
addpath(cd_org)

path_tools = '/usr/local/MATLAB_Tools';
set_ft(path_tools)

%%
% cd('/MEG_data/epilepsy')
% [analyses_dir] = uigetdir; % choose, /MEG_data/epilepsy/xxx/211015/analyses
% cd(analyses_dir)

%% MR_org
% cd('/MEG_data/MRI_database/epilepsy/WAGNER_Sarah_sagT1/DICOM')
cd('/MEG_data/MRI_database/epilepsy/')
set_spm(path_tools)
dicomfile_MR = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

set_ft(path_tools)
dicom_MR = ft_read_mri(dicomfile_MR);
ft_sourceplot([], dicom_MR);

%% Saving data as nii
tkz = tokenize(dicomfile_MR,'/');
disp(tkz');

% disp('/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare/Export/Parra_J')
savedir = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare/Export/';
% savedir = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom_prepare/Export/';
% savedir = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare/Export/Parra_J';
cd(savedir)
subname = input('enter subject name last_first:','s');
if exist(fullfile(savedir,subname), 'file') == 0, mkdir(fullfile(savedir,subname)); end
% T1 = ft_read_mri('T1_BS.nii');

subdir = fullfile(savedir,subname);

%%
cd('/MEG_data/epilepsy/')
disp('select: /xxxx/211015/report/Spikes')
reportdir = uigetdir; 

% reportdir = '/MEG_data/epilepsy/parra_jocelyn/211015/report/Spikes';
cd(reportdir);

%%
T1_nii = uigetfile ({'*.nii','T1 (*.nii)'},'select T1 nii');
T1 = ft_read_mri(T1_nii); ft_sourceplot([], T1);
dip_nii = uigetfile ({'*.nii','dip (*.nii)'},'select dipole nii');
dip = ft_read_mri(dip_nii); ft_sourceplot([], dip);

%%
cd(subdir)
nii_name = 'dip.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, dip);
ft_sourceplot([], dip);

%%
cd(subdir)
nii_name = 'dicom_MR.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, dicom_MR);

nii_name = 'T1.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, T1);
ft_sourceplot([], T1);

%%
% T1_flirt = ft_read_mri('T1_flirt.nii');

%% SPM coreg, estimate and reslice
cd(subdir)
addpath(cd_org)
set_spm(path_tools)
% set_spm

nii_filename1 = 'dicom_MR.nii'; % 256x256x150
nii_filename2 = 'T1.nii';       % 256x256x256
nii_filename3 = 'dip.nii';      % 256x256x256

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[nii_filename1,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[nii_filename2,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[nii_filename3,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch);

%%
close all
set_ft(path_tools)
rT1 = ft_read_mri('rT1.nii');
r_dip = ft_read_mri('rdip.nii');
ft_sourceplot([], r_dip);
ft_sourceplot([], rT1);

%% Export overlayed dipole and T1
ovT1_dip = rT1;
tmp = r_dip.anatomy+rT1.anatomy;
ovT1_dip.anatomy = tmp;

nii_name = 'T1_dip_overlay.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, ovT1_dip);
ft_sourceplot([], ovT1_dip);

%%
[a,b] = fileparts(dicomfile_MR);
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/FT_fucntions/External/Miscellaneous');
% d = rdir(fullfile(a,'EXP*'));
d = rdir(fullfile(a,'dicom-*'));

metadata_all = [];
for i=1:length(d)
    metadata = dicominfo(d(i).name);
%     disp(d(i).name)
    metadata_all{i}=metadata;
%     metadata.SliceLocation
end

%%
outd = fullfile(subdir,'dicom','T1'); if exist(outd, 'file') == 0, mkdir(outd); end

cd(subdir)
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
    
    metadata  = metadata_all{i};    
    dicomImage = dicomread(fullfile(outd,dicome_name));
    dicomInfo = dicominfo(fullfile(outd,dicome_name));
    dicomInfo.SliceLocation = metadata.SliceLocation;
    dicomInfo.ImagePositionPatient = metadata.ImagePositionPatient;
    dicomInfo.ImageOrientationPatient = metadata.ImageOrientationPatient;
    dicomInfo.PixelSpacing = metadata.PixelSpacing;
    dicomwrite(dicomImage, fullfile(outd,dicome_name), dicomInfo, 'CreateMode', 'copy');

end

outd = fullfile(subdir,'dicom','dip'); if exist(outd, 'file') == 0, mkdir(outd); end
for i=1:size(r_dip.anatomy,3)
    I = int16(squeeze(r_dip.anatomy(:,:,i)));
    if i < 11
        dicome_name = ['EXP000',num2str(i-1)];
    elseif i < 101
        dicome_name = ['EXP00',num2str(i-1)];
    else
        dicome_name = ['EXP0',num2str(i-1)];
    end
    %     dicomwrite(I,dicome_name, metadata_all{i});
    dicomwrite(I,fullfile(outd,dicome_name),metadata_all{i});
    
    metadata  = metadata_all{i};
    dicomImage = dicomread(fullfile(outd,dicome_name));
    dicomInfo = dicominfo(fullfile(outd,dicome_name));
    dicomInfo.SliceLocation = metadata.SliceLocation;
    dicomInfo.ImagePositionPatient = metadata.ImagePositionPatient;
    dicomInfo.ImageOrientationPatient = metadata.ImageOrientationPatient;
    dicomInfo.PixelSpacing = metadata.PixelSpacing;
    dicomwrite(dicomImage, fullfile(outd,dicome_name), dicomInfo, 'CreateMode', 'copy');
end


% Overlayed T1 and Dipople
outd = fullfile(subdir,'dicom','T1_dip'); if exist(outd, 'file') == 0, mkdir(outd); end

cd(subdir)
for i=1:size(ovT1_dip.anatomy,3)
    I = int16(squeeze(ovT1_dip.anatomy(:,:,i)));
    if i < 11
        dicome_name = ['EXP000',num2str(i-1)];
    elseif i < 101
        dicome_name = ['EXP00',num2str(i-1)];
    else
        dicome_name = ['EXP0',num2str(i-1)];
    end
    %     dicomwrite(I,dicome_name, metadata_all{i});
    dicomwrite(I,fullfile(outd,dicome_name), metadata_all{i});
    
    metadata  = metadata_all{i};
    dicomImage = dicomread(fullfile(outd,dicome_name));
    dicomInfo = dicominfo(fullfile(outd,dicome_name));
    dicomInfo.SliceLocation = metadata.SliceLocation;
    dicomInfo.ImagePositionPatient = metadata.ImagePositionPatient;
    dicomInfo.ImageOrientationPatient = metadata.ImageOrientationPatient;
    dicomInfo.PixelSpacing = metadata.PixelSpacing;
    dicomwrite(dicomImage, fullfile(outd,dicome_name), dicomInfo, 'CreateMode', 'copy');
end

%% CHECKING DICOM outputs
set_spm(path_tools)
dicomfile_T1_check = spm_select(1,'.*','Select T1 dicome file, e.g. EXP0000');

set_ft(path_tools)
T1_check = ft_read_mri(dicomfile_T1_check);
ft_sourceplot([], T1_check);

nii_name = 'dicom_out_T1.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, T1_check);

close all
set_spm(path_tools)
dicomfile_dip_check = spm_select(1,'.*','Select dipole dicome file, e.g. EXP0000');

set_ft(path_tools)
dip_check = ft_read_mri(dicomfile_dip_check);
ft_sourceplot([], dip_check);

nii_name = 'dicom_out_dip.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, dip_check);


close all
set_spm(path_tools)
dicomfile_T1_Ldip_check = spm_select(1,'.*','Select T1_dipole dicome file, e.g. EXP0000');

set_ft(path_tools)
T1_dip_check = ft_read_mri(dicomfile_T1_Ldip_check);
ft_sourceplot([], T1_dip_check);

nii_name = 'dicom_out_T1_dip.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, T1_dip_check);

%%
set_spm(path_tools)
spm_check_registration

%%
function set_ft(path_tools)

restoredefaultpath
% ft_path = fullfile(path_tools,'/fieldtrip-20161201');
ft_path = fullfile(path_tools,'/fieldtrip_20190419');
addpath(fullfile(path_tools,'/fieldtrip_20190419/external/spm8'));
addpath(ft_path);
ft_defaults

end

function set_spm(path_tools)

restoredefaultpath
spm_path = fullfile(path_tools,'/spm12');
addpath(genpath(spm_path));
spm_get_defaults

end
