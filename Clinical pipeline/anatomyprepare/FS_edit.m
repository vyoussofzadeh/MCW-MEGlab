clear, clc, close all,

% Author, Vahab Youssof Zadeh, 2021-22
% updates: 7/26/22

%%
cd_org = '/MEG_data/MCW_pipeline/Anatomyprepare';
addpath(cd_org)
path_tools = '/usr/local/MATLAB_Tools';

%%
disp('Edit has to be done at the FS folder dir!')
% disp('1: Edit from MRI_database')
disp('2: Edit from FS subject dir')
% in_sel = input(':');
in_sel = 2;

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
% command = ['tkmedit ', name, ' T1.mgz -surfs'];
% 
% disp('1: Yes');
% disp('2: No');
% fschkask = input('Edit FS?');
% if fschkask ==1
%     clc, close all,
%     disp('run this in command ...'),
%     disp(['cd ', pwd]);
%     disp(command)
% end

%%
clc
disp('source:')
disp('https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/ControlPoints_freeview')
disp('========');

disp(['freeview -v ', fullfile(pathstr, name), '/mri/brainmask.mgz \'])
disp([' -f ', fullfile(pathstr, name),'/surf/lh.white:edgecolor=blue \'])
disp([fullfile(pathstr, name),'/surf/lh.pial:edgecolor=red \'])
disp([fullfile(pathstr, name),'/surf/rh.white:edgecolor=blue \'])
disp([fullfile(pathstr, name),'/surf/rh.pial:edgecolor=red']);

%%
disp('============')
disp('File/newsetpoint: control.dat')
disp('Add, click (and remove, shift_click) ctrl points')
disp('save control ppoints in the  FS tmp folder')
disp('============')
disp('then try ...')
disp(['cd ', fullfile(pathstr)]);
disp(['recon-all -autorecon2-cp -autorecon3 -subjid ', fullfile(pathstr,name)])