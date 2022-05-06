restoredefaultpath
addpath(cd_org)
spm_path = fullfile(path_tools,'/spm12');
addpath(genpath(spm_path));
spm_get_defaults