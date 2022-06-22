%% The bgros MEG pilot data

% MEG analysis pipeline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 06/10/2022

clear; clc, close('all'); warning off

%% Analysis flags
flag = [];
flag.plottriggers = 1;
flag.anatomy_check = 1;
flag.dicsanalysis = 0;
flag.connanalysis = 1;

%% Paths
addpath('./run')
Run_setpath

%% Anatomy
subj = 'pilot1';
run = 1;
datafilename = ['Task_run', num2str(run), '_raw.fif'];
datafile = fullfile(indir, datafilename);
subdir = fullfile(outdir, subj,'preprocess'); % output dir
load(fullfile(subdir,['data_run', num2str(run), '_cond',num2str(1), '.mat']))

mridir = '/group/bgross/work/CIDMEG/analysis/anatomy/Sub001/T1';
mrifile = 'sub-CIDFMRI001_ses-1_acq-mprage_T1w.nii.gz';
Run_anatomy

%%
close all
savedir = '/group/bgross/work/CIDMEG/analysis/process/pilot1/process';
ddd = rdir(fullfile(savedir, '/**/wPLI*.mat'));

clear ID conn_val
for i=1:length(ddd)
    tkz = tokenize(ddd(i).name,'/');
    conn = load(ddd(i).name);
    conn_val(:,:,:,i) = conn.source_conn1.wpli_debiasedspctrm;
    ID{i} = tkz{end}(1:end-4);
end
mconn_val = nanmean(conn_val,4);


source_conn1 = conn.source_conn1;
source_conn1.wpli_debiasedspctrm = mconn_val;
par = conn.par;
vs_roi1 = conn.vs_roi1;

cd('/group/bgross/work/CIDMEG/analysis/process/pilot1/process/group')
Run_networkanalysis %- network analysis (network degrees)
Run_vis_netconn_analysis %- mapping
