clear, clc, close all,

% Author, Vahab Youssof Zadeh, 2021-22
% updates: 12/16/21, 7/26/22

%%
cd_org = '/MEG_data/MCW_pipeline/Anatomyprepare';
addpath(cd_org)
path_tools = '/usr/local/MATLAB_Tools';

%%
disp('1: Check from MRI_database')
disp('2: Check from FS subject dir')
in_sel = input(':');

switch in_sel
    case 1
        indir = '/MEG_data/MRI_database/epilepsy/';
    case 2
        indir = '/usr/local/freesurfer6/subjects';
end

cd(indir)
FSdir = uigetdir;
[pathstr, name] = fileparts(FSdir);
cd(pathstr)

%%
command = ['tkmedit ', name, ' T1.mgz -surfs'];

disp('1: Yes');
disp('2: No');
fschkask = input('Check FS?');
if fschkask ==1
    clc, close all,
    disp('run this in command ...'),
    disp(['cd ', pwd]);
    disp(command)
end

% /MEG_data/MRI_database/epilepsy/RAPEY_Ward_Jacob/FSrecon_110921
% cd /MEG_data/MRI_database/epilepsy/RAPEY_Ward_Jacob/
% tkmedit FSrecon_110921 T1.mgz -surfs

%%


