clc;
clear;
close all;

% runMaxFilter (t)SSS for MEGIN data
% Writtern by MCW MEG group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>

% Update, SSS analysis was added for emptyroom data, 03/24/25
% Update: Pipeline was created 09/11/24

%% Remove dynamically added Java class paths
dynamicPaths = javaclasspath('-dynamic');
for idx = length(dynamicPaths):-1:1
    javarmpath(dynamicPaths{idx});
end
clear java;

%%
do_analysis = [];
do_analysis.ecg = 1;
do_analysis.eog = 0;
do_analysis.ongoing = 1;

%% Add necessary paths
pathList = {'/MEG_data/megclinic',...
    '/usr/local/MATLAB_Tools/brainstorm3',...
    '/usr/local/MNE-2.7.0-3106-Linux-x86_64/share/matlab', ...
    '/usr/local/MATLAB_Tools/mne', ...
    '/neuro/bin',...
    '/MEG_data/MCW_pipeline/Preprocess/func'};
for p = pathList
    addpath(genpath(char(p)));
end

subjIDpath = '/MEG_data/LAB_MEMBERS/Vahab/savingdatatest/2';
outDirpath = '/MEG_data/LAB_MEMBERS/Vahab/savingdatatest/2';
cleanup_MaxFilter_tSSS_forWrapper(subjIDpath, outDirpath)
