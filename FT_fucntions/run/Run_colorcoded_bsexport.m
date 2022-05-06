Scouts = atlas.Scouts;
TessNbVertices = atlas.TessNbVertices;
nScouts = length(Scouts);
colr = hsv(nScouts);
colr = hot(nScouts);
% colr = viridis(nScouts);
addpath('/data/MEG/Vahab/Github/MCW-MEGlab/tools/helpful_tools');
colr = viridis(nScouts);

% left and right regions
% idx_left = 1:2:nScouts;
% idx_right = 2:2:nScouts;
% idx_all = [idx_left, idx_right];

% thre = input('enter threshold value:');
thre = 0;

read_atlas_lang = atlas;

[l, idx] = sort((datain),'descend');
datain(isnan(datain))=0;
roiid(idx(1:nroi))'

k=1;
read_atlas_lang1 = [];
for i=1:nScouts
    if abs(datain(i)) >= thre*max(datain)
        read_atlas_lang1.Scouts(k) = read_atlas_lang.Scouts(idx(i));
        read_atlas_lang1.Scouts(k).Color = colr(end-k+1,:);
        k=k+1;
        %                 disp(l(i))
    else
        read_atlas_lang1.Scouts(i) = read_atlas_lang.Scouts(idx(i));
        read_atlas_lang1.Scouts(idx(i)).Color = colr(1,:);
        disp(i)
    end
end
read_atlas_lang1.Name = savetag;
read_atlas_lang1.TessNbVertices = read_atlas_lang.TessNbVertices;
save(fullfile(savetag),'-struct', 'read_atlas_lang1'),
