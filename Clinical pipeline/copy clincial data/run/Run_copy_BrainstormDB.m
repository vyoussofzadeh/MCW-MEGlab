%% Copying brainstorm files
clc
idx = [];
idx.matched = data_idx(ib);

clc
clear datafile_sub
k=1;
for i=1:length(name_list.newmatched)
    fname = name_list.newmatched{i};
    idx.matched2 = find(contains(name_list.all,fname)==1);
    if ~isempty(idx.matched2)
        savepfolder = fullfile(savedir, name_list.all{idx.matched2}, 'brainstorm_db');
        if exist(savepfolder, 'file') == 0, mkdir(savepfolder), end
        cd(savepfolder)
        d = dir(fullfile(indir,fname, '/brainstorm_db'));
        if ~isempty(d) && ~exist(fullfile(savepfolder,'data'),'file') && ~exist(fullfile(savepfolder,'anat'),'file')
            disp([num2str(i),': copying: ', fname])
            disp(fullfile(indir,fname, '/brainstorm_db/data')), disp('to')
            disp(fullfile(savepfolder, 'data'))
            copyfile(fullfile(indir,fname, '/brainstorm_db/data'), fullfile(savepfolder, 'data'))
            copyfile(fullfile(indir,fname, '/brainstorm_db/anat'), fullfile(savepfolder, 'anat'))
        else
            disp([num2str(i), ': check data: ', fname])
        end
    end
end