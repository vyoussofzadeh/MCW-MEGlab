
clc, clear, close all,

%%
% cd_org = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare';
path_tools = '/usr/local/MATLAB_Tools';
set_ft(path_tools)

%%
% cd('/MEG_data/epilepsy')
% [analyses_dir] = uigetdir; % choose, /MEG_data/epilepsy/xxx/211015/analyses
% cd(analyses_dir)

%% MR_org
% cd('/MEG_data/MRI_database/epilepsy/WAGNER_Sarah_sagT1/DICOM')
cd('/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/dicom prepare/CD')
set_spm(path_tools)
dicomfile_MR = spm_select(1,'.*','Select one dicome file, e.g. EXP0000');

set_ft(path_tools)
dicom_MR = ft_read_mri(dicomfile_MR);
ft_sourceplot([], dicom_MR);

%%
nii_name = 'CHECK.nii';
cfg                 = [];
cfg.filename        = nii_name;
cfg.filetype        = 'nifti';
cfg.parameter       = 'anatomy';
ft_volumewrite(cfg, dicom_MR);
ft_sourceplot([], dicom_MR);


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
