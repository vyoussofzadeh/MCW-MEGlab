function set_spm(path_tools)

restoredefaultpath
spm_path = fullfile(path_tools,'spm12');
addpath(genpath(spm_path));
spm_get_defaults

end