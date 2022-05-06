coor_file = '/data/MRI/ECP/resting/scripts/centers.tsv';

s = tdfread(coor_file);
mni_atlas = [s.x,s.y,s.z];
conn = rand(length(mni_atlas),length(mni_atlas));

idx = randperm(length(mni_atlas));
idx_sel = idx(1:10); % 10 random connection nodes

conn_mesh_display('', '', '', mni_atlas(idx_sel,:), tedge(idx_sel,idx_sel), .2);
