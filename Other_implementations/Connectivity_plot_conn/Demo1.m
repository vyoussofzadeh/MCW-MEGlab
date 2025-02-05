
conn_path = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Conn/conn';
addpath(conn_path);

spmpath = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/SPM/spm12_2021/spm12';
addpath(spmpath)

coor_file = 'centers.tsv';

s = tdfread(coor_file);
mni_atlas = [s.x,s.y,s.z];
conn = rand(length(mni_atlas),length(mni_atlas));

idx = randperm(length(mni_atlas));
idx_sel = idx(1:10); % 10 random connection nodes

conn_mesh_display('', '', '', mni_atlas(idx_sel,:), conn(idx_sel,idx_sel), .2);
