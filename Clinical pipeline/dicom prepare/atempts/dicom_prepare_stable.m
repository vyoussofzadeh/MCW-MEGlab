clc, clear, close all,

%%
% cd_org = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare';
path_tools = '/usr/local/MATLAB_Tools';
set_ft(path_tools)

%%  Saving data as nii
savedir = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare/SWagner';
cd(savedir)
% T1 = ft_read_mri('T1_BS.nii');

%%
reportdir = '/MEG_data/epilepsy/wagner_sarah/210621/report/Spikes';
cd(reportdir);
T1 = ft_read_mri('T1.nii');
L_dip = ft_read_mri('run2-3-4-5_dipoles_90gof_Lsensors.nii');

%%
cd(savedir)
nii_name = 'L_dip.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, L_dip);
ft_sourceplot([], L_dip);

%% MR_org
cd('/MEG_data/MRI_database/epilepsy/WAGNER_Sarah_sagT1/DICOM')

set_spm(path_tools)
dicomfile_MR = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

set_ft(path_tools)
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
cd(savedir)

set_spm(path_tools)
nii_filename1 = 'dicom_MR.nii'; % 256x256x150
nii_filename2 = 'T1.nii';       % 256x256x256
nii_filename3 = 'L_dip.nii';       % 256x256x256


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
rL_dip = ft_read_mri('rL_dip.nii');
ft_sourceplot([], rL_dip);
ft_sourceplot([], rT1);

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
set_spm(path_tools)
dicomfile_T1_check = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

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
dicomfile_Ldip_check = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

set_ft(path_tools)
ldip_check = ft_read_mri(dicomfile_Ldip_check);
ft_sourceplot([], ldip_check);

nii_name = 'dicom_out_ldip.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, ldip_check);

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
