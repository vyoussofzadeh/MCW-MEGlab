% clear; clc, close('all'); warning off
% 
% bs_path = '/opt/matlab_toolboxes/brainstorm3';
% set(0,'DefaultFigureWindowStyle' , 'normal')
% addpath(bs_path);
% brainstorm
% addpath(genpath('./functions'));
% 
% %%
% % datafile = '/data/MEG/Vahab/test_data/raghavan_manoj2/brainstorm_db/anat/raghavan_manoj2';
% % datafile = '/data/MEG/Vahab/Scripts/Vahab/Test_MEG2/bednar_peggy/brainstorm_db/anat/bednar_peggy';
% % datafile = '/data/MEG/Clinical/MEG/dougherty_danielle/brainstorm_db/anat/dougherty_danielle';
% % datafile = '/data/MEG/Vahab/test_data/busby_michael/brainstorm_db/anat/busby_michael';
% % datafile = '/data/MEG/Vahab/test_data/krier_henry/brainstorm_db/anat/krier_henry';
% % datafile = '/data/MEG/Vahab/test_data/anderson_mark/brainstorm_db/anat/anderson_mark';
% % datafile = '/data/MEG/Vahab/test_data/brainstorm_db/anat/meier_kayla';
% datafile = '/data/MEG/Vahab/test_data/berg_gwyneth/brainstorm_db/anat/berg_gwyneth';
% 
% 
% disp(datafile)
% d = rdir(fullfile(datafile,'subjectimage*.mat'));
% % d = rdir(fullfile(datafile,'*ras.mat'));
% if ~isempty(d)
%     sMRI = d.name;
%     cd(datafile)
%     cd ..
%     OutputFile = fullfile(pwd,'T1.nii');
%     out_mri_nii(sMRI, OutputFile);
% end

% clear; clc, close('all'); warning off

%%
% bs_path = '/opt/matlab_toolboxes/brainstorm3';
% addpath(bs_path);
% brainstorm
% addpath(genpath('./functions'));

close all
set(0,'DefaultFigureWindowStyle','normal')
addpath('/opt/matlab_toolboxes/brainstorm3')
brainstorm

[a,b] = fileparts(cfg_main.mripfile);

anatfile = fullfile(cfg_main.mripfile);

%%
% datafile = '/data/MEG/Vahab/test_data/raghavan_manoj2/brainstorm_db/anat/raghavan_manoj2';
% datafile = '/data/MEG/Vahab/Scripts/Vahab/Test_MEG2/bednar_peggy/brainstorm_db/anat/bednar_peggy';
% datafile = '/data/MEG/Clinical/MEG/dougherty_danielle/brainstorm_db/anat/dougherty_danielle';
% datafile = '/data/MEG/Vahab/test_data/busby_michael/brainstorm_db/anat/busby_michael';
% datafile = '/data/MEG/Vahab/test_data/krier_henry/brainstorm_db/anat/krier_henry';
% datafile = '/data/MEG/Vahab/test_data/anderson_mark/brainstorm_db/anat/anderson_mark';


disp(anatfile)
d = rdir(fullfile(a,'subjectimage*.mat'));
% d = rdir(fullfile(datafile,'*ras.mat'));
if ~isempty(d)
    sMRI = d.name;
    cd(a)
%     cd ..
    OutputFile = fullfile(pwd,'T1.nii');
    out_mri_nii(sMRI, OutputFile);
end

%%
brainstorm stop
rmpath(genpath('/opt/matlab_toolboxes/brainstorm3'));
