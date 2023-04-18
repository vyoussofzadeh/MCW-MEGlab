restoredefaultpath
addpath(cd_org)
% spm_path = fullfile(path_tools,'/spm12');
% addpath(genpath(spm_path));
% spm_get_defaults
% 
% 
% addpath('/usr/local/MATLAB_Tools/spm12')
% spm_defaults

spmpath = '/usr/local/MATLAB_Tools/spm12';
addpath(spmpath)
% spm_get_defaults
spm('defaults', 'PET');