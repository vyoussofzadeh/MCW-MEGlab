clear; clc, close('all'); warning off

%%
%- Input dir
fileDir = '/MEG_data/Vahab/Github/MCW-MEGlab/FT/Clinical pipeline/sessionlog';
% fileDir = 'D:\excel_folder';
outfile = 'sessionlog/master_session.xlsx';
delete(outfile)
addpath(fileDir);
fileNames = dir(fileDir);
fileNames = {fileNames.name};
fileNames = fileNames(cellfun(@(f)contains(f,'.xlsx'),fileNames));
% for f = 1:numel(fileNames)
%     fTable = readtable(fileNames{f});
%     writetable(fTable,outfile,'Sheet',fileNames{f});
% end

T = [];
for f = 1:numel(fileNames)
    fTable = readtable(fileNames{f},'sheet','Sheet1');
    T = [T;fTable];
end
writetable(T,outfile,'sheet','patientdata');
disp('completed!')