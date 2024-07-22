% restoredefaultpath
% addpath(cd_org)
% spm_path = fullfile(path_tools,'/spm12');
% addpath(genpath(spm_path));
% spm_get_defaults

restoredefaultpath
addpath(cd_org)
% spm_path = ('/data/MEG/Vahab/Github/MCW_MEGlab/tools/SPM');
% addpath(genpath(spm_path));
% spm_get_defaults

spmpath = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/SPM/spm12_2021/spm12';
addpath(spmpath)
% spm_get_defaults
spm('defaults', 'PET');

% % spmpath = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/SPM/spm12_2021/spm12';
% addpath(spmpath)
% % spm_get_defaults
% spm('defaults', 'PET');