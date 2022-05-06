clear
addpath('/data/MEG/Vahab/Github/MCW-MEGlab/tools/Conn/conn');
h = get(0, 'Children');
if isempty(findobj(h,'tag','CONN functional connectivity toolbox'))
    conn;
end
path_tools = '/data/MEG/Vahab/Github/MCW-MEGlab/tools';
spm_path = fullfile(path_tools,'SPM/spm12');
addpath(genpath(spm_path))
spm_get_defaults

addpath('/data/MEG/Vahab/Github/MCW-MEGlab/FT/functions/External/brewermap')

tedge= 1*rand(3,3);
mni = [1,20,10; -10,1,-20; 60,2,25];
conn_mesh_display('', '', '', mni, tedge, .2);

%%
coor_file = '/data/MRI/ECP/resting/scripts/centers.tsv';

s = tdfread(coor_file);
mni_atlas = [s.x,s.y,s.z];
conn = rand(length(mni_atlas),length(mni_atlas));

idx = randperm(length(mni_atlas));
idx_sel = idx(1:3);

conn_mesh_display('', '', '', mni_atlas(idx_sel,:), tedge, .2);
%%
