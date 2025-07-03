%% ECP Semantic Decision Task Dataset, Medical College of Wisconsin

% Script: BS Process (Laterality Analysis)
% Project: ECP_SD
% Written by: Vahab Youssof Zadeh

clear; clc; close all; warning off;

%% Paths
restoredefaultpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/run')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/function')
Run_setpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other/')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/tools/helpful_tools/daviolinplot/daboxplot');

%% 72 patients
sub_MF_pt = {'EC1007';'EC1014';'EC1020';'EC1028';'EC1029';'EC1030';'EC1042';'EC1047';'EC1048';'EC1060';'EC1067';'EC1068';'EC1072';'EC1073';'EC1074';'EC1075';'EC1077';'EC1078';'EC1079';'EC1081';'EC1082';'EC1097';'EC1099';'EC1100';'EC1101';'EC1102';'EC1103';'EC1106';'EC1109';'EC1113';'EC1114';'EC1115';'EC1116';'EC1117';'EC1118';'EC1121';'EC1122';'EC1123';'EC1125';'EC1129';'EC1130';'EC1133';'EC1134';'EC1135';'EC1138';'EC1139';'EC1141';'EC1142';'EC1144';'EC1150';'EC1153';'EC1154';'EC1156';'EC1157';'EC1158';'EC1159';'EC1161';'EC1162';'EC1165';'EC1167';'EC2029';'EC2038';'EC2045';'EC2047';'EC2052';'EC2054';'EC2072';'EC2074';'EC2083';'EC2085';'EC2090';'EC2109';'EC2112';'EC2114'};

%% Demographic Details
patn_MEGfMRI_neuropsych = ecpfunc_sum_neuropsych(sub_MF_pt);

%%
outdir = '/data/MEG/Research/ECP/Behavioural/update_vzadeh_MEGlist_070225';
fname = fullfile(outdir,'patn_MEGfMRI_neuropsych.xlsx');
writetable(patn_MEGfMRI_neuropsych, fname);   % keep it simple
