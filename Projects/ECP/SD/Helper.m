clc
cd('/data/MEG/Research/ECP/Semantic_Decision/BS_database/data/Group_analysis')
d = rdir('./ec*/results_abs_*.mat');

for j = 1:length(d)
    disp(d(j).name)
%     pause,
    delete(d(j).name)
end
