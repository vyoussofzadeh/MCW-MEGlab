addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/brewermap')

%% SPM
spmpath = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/SPM/spm12_2021/spm12';
addpath(spmpath)
spm_get_defaults

%% Conn
conn_path = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Conn/conn';
addpath(genpath(conn_path));

%%
cd('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/neurovault_MMP/MMP 1.0 MNI projections')

% HCP = ft_read_mri('MMP_in_MNI_symmetrical.nii');
% unique(HCP.anatomy)


conn_mesh_display('MMP_in_MNI_symmetrical.nii');

%%
% cd('/data/MEG/Vahab/Github/surface-atlas-FieldTrip-source-analysis/HCP')
% HCP = ft_read_mri('MMP_in_MNI_corr.nii');
% unique(HCP.anatomy)
% conn_mesh_display('MMP_in_MNI_corr.nii');