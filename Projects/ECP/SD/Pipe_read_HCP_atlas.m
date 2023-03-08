%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (HCP atlas)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 01/11/2023

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
cd('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD')
addpath('./run')
Run_setpath
addpath('./data')
addpath('./run')
% addpath(genpath('./functions'))
addpath(genpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions'))
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/Colormaps-from-MatPlotLib2.0')

results_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results';
addpath(ft_path); ft_defaults

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new')

s_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/glasser atlas';

%% Read manual networks for LI analysis (Networks were provided by Joe)
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["parcel_num", "num_vertices", "roi_num", "roi_name", "roi_count", "roi_percentage", "bilateral_parcel_name", "bilateral_roi_name", "bilateral_roi_count", "bilateral_num_vertices", "bilateral_roi_count_max", "bilateral_roi_percentage", "manual_roi_name"];
opts.VariableTypes = ["double", "double", "double", "string", "double", "double", "string", "string", "double", "double", "double", "double", "categorical"];
opts = setvaropts(opts, 7, "WhitespaceRule", "preserve");
opts = setvaropts(opts, [4, 7, 8, 13], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
manualglasserROISSDTD2020 = readtable("/group/jbinder/ECP/MEG/laterality_index/manual_glasser_ROIS_SDTD_2020.tsv", opts);

clear opts
glass_net_regions = unique(manualglasserROISSDTD2020.roi_name);

glass_net_L_idx = [];  glass_net_R_idx = [];
glass_net_L_label = [];  glass_net_R_label = [];
k1 = 1; k2 = 1;
for i=1:length(glass_net_regions)
    disp(glass_net_regions(i))
    if contains(glass_net_regions(i),'L')
        %         glass_roi_L_idx = [glass_roi_L_idx, i];
        idx = find(contains(manualglasserROISSDTD2020.roi_name, glass_net_regions(i))==1);
        glass_net_L_idx{k1} = idx;
        glass_net_L_label{k1} = manualglasserROISSDTD2020.bilateral_parcel_name(idx);
        
        tmp  = glass_net_L_label{k1};
        for k=1:length(tmp), tmp{k} = ['L_', tmp{k}]; end
        glass_net_L_label{k1} = tmp;
        
        k1 = k1 + 1;
        disp(manualglasserROISSDTD2020.roi_name(idx))
    elseif contains(glass_net_regions(i),'R')
        %         glass_roi_R_idx = [glass_roi_L_idx, i];
        idx = find(contains(manualglasserROISSDTD2020.roi_name, glass_net_regions(i))==1);
        glass_net_R_idx{k2} = idx;
        glass_net_R_label{k2} = manualglasserROISSDTD2020.bilateral_parcel_name(idx);
        tmp  = glass_net_R_label{k2};
        for k=1:length(tmp), tmp{k} = ['R_', tmp{k}]; end
        glass_net_R_label{k2} = tmp;
        
        k2 = k2 + 1;
        disp(manualglasserROISSDTD2020.roi_name(idx))
    end
    %     pause
end

save(fullfile(s_dir, 'LI_glasser_manual_net_12.mat'),  'glass_net_regions', 'glass_net_L_idx','glass_net_R_idx', 'glass_net_R_label', 'glass_net_L_label');

cd(s_dir)
