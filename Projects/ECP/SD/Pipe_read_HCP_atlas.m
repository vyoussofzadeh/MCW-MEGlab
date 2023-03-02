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

%%
% atlas = load('/data/MEG/Vahab/Github/MCW_MEGlab/FT/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_corr.nii_362_updated.mat');
% 
% clear rois region region_lr
% for i=1:length(atlas.Scouts)
%     rois{i} = atlas.Scouts(i).Label;
%     region{i} = atlas.Scouts(i).Region;
%     region_lr{i} = atlas.Scouts(i).Region(1);
% end
% idx_L = []; idx_R = []; idx_clr_L = []; idx_clr_R = [];
% for i = 1:length(rois)
%     idx = [];
%     if find(strfind(rois{i}, 'L_')==1)
%         idx_L = [idx_L; i];
%     elseif find(strfind(rois{i}, 'R_')==1)
%         idx_R = [idx_R; i];
%     end
% end
% length(idx_R)
% length(idx_L)
% 
% [groups_labels, groups] = do_read_HCP_labels_bs;

%% Plot HCP atlas
% cfg = [];
% cfg.atlas = atlas;
% cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
% cfg.sel = 'whole'; % 'whole', 'left', 'right', 'roi';
% cfg.lat_index = [idx_L, idx_R];
% cfg.rois = rois;
% do_plot_HCP_atlas(cfg);

%% Bilateral (symmetrical) ROIs for LI analysis (ROIs were provided by Joe)
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["parcel_num", "parcel_name", "num_vertices", "roi_num", "roi_name", "roi_count", "roi_percentage", "bilateral_parcel_name", "bilateral_roi_name", "bilateral_roi_count", "bilateral_num_vertices", "bilateral_roi_count_max", "bilateral_roi_percentage"];
opts.VariableTypes = ["double", "string", "double", "double", "string", "double", "double", "string", "categorical", "double", "double", "double", "double"];
opts = setvaropts(opts, [2, 8], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 5, 8, 9], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
bilateralglasserlateral = readtable("/group/jbinder/ECP/MEG/laterality_index/bilateral_glasser_lateral.tsv", opts);
clear opts

glass_roi_L_idx = []; 
glass_roi_R_idx = [];
glass_roi_R_name = [];
glass_roi_L_name = [];
for i=1:length(bilateralglasserlateral.roi_name)
    if contains(bilateralglasserlateral.roi_name(i),'L_')
        glass_roi_L_idx = [glass_roi_L_idx, i]; 
        glass_roi_L_name = [glass_roi_L_name; bilateralglasserlateral.parcel_name(i)];
    elseif contains(bilateralglasserlateral.roi_name(i),'R_')
        glass_roi_R_idx = [glass_roi_R_idx, i];
        glass_roi_R_name = [glass_roi_R_name; bilateralglasserlateral.parcel_name(i)];
    end
end

s_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/glasser atlas';
glass_roi_name = [glass_roi_L_name; glass_roi_R_name];

%- Lookup
glass_roi_L_idx = []; 
glass_roi_R_idx = [];
for i = 1:length(glass_roi_L_name)
    if ~isempty(find(contains(bilateralglasserlateral.parcel_name, glass_roi_L_name(i))==1, 1))
        idx = find(contains(bilateralglasserlateral.parcel_name, glass_roi_L_name(i))==1);
        glass_roi_L_idx = [glass_roi_L_idx, idx];
        disp([bilateralglasserlateral.parcel_name(idx), glass_roi_L_name(i)])
    end
end
for i = 1:length(glass_roi_R_name)
    if ~isempty(find(contains(bilateralglasserlateral.parcel_name, glass_roi_R_name(i))==1, 1))
        idx = find(contains(bilateralglasserlateral.parcel_name, glass_roi_R_name(i))==1);
        glass_roi_R_idx = [glass_roi_R_idx, idx];
        disp([bilateralglasserlateral.parcel_name(idx), glass_roi_R_name(i)])
    end
end
save(fullfile(s_dir, 'LI_glasser_bilateral_rois_90.mat'), 'glass_roi_name','glass_roi_R_idx','glass_roi_L_idx');


%% Read bilateral networks for LI analysis (Networks were provided by Joe)
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["parcel_num", "parcel_name", "num_vertices", "roi_num", "roi_name", "roi_count", "roi_percentage", "bilateral_parcel_name", "bilateral_roi_name", "bilateral_roi_count", "bilateral_num_vertices", "bilateral_roi_count_max", "bilateral_roi_percentage"];
opts.VariableTypes = ["double", "string", "double", "double", "string", "double", "double", "string", "categorical", "double", "double", "double", "double"];
opts = setvaropts(opts, [2, 8], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 5, 8, 9], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
% bilateralglasserROISSDTD2020 = readtable("/group/jbinder/ECP/MEG/laterality_index/bilateral_glasser_ROIS_SDTD_2020.tsv", opts);
bilateralglasserROISSDTD2020 = readtable("/group/jbinder/ECP/MEG/laterality_index/manual_glasser_ROIS_SDTD_2020.tsv", opts);

clear opts

glass_net_regions = unique(bilateralglasserROISSDTD2020.roi_name);


glass_net_L_idx = [];  glass_net_R_idx = [];
glass_net_L_label = [];  glass_net_R_label = [];
k1 = 1; k2 = 1;
for i=1:length(glass_net_regions)
    disp(glass_net_regions(i))
    if contains(glass_net_regions(i),'L')
%         glass_roi_L_idx = [glass_roi_L_idx, i];
        idx = find(contains(bilateralglasserROISSDTD2020.roi_name, glass_net_regions(i))==1);
        glass_net_L_idx{k1} = idx; 
        glass_net_L_label{k1} = bilateralglasserROISSDTD2020.parcel_name(idx);
        k1 = k1 + 1;
        disp(bilateralglasserROISSDTD2020.roi_name(idx))
    elseif contains(glass_net_regions(i),'R')
%         glass_roi_R_idx = [glass_roi_L_idx, i];
        idx = find(contains(bilateralglasserROISSDTD2020.roi_name, glass_net_regions(i))==1);
        glass_net_R_idx{k2} = idx;
        glass_net_R_label{k2} = bilateralglasserROISSDTD2020.parcel_name(idx);
        k2 = k2 + 1;
        disp(bilateralglasserROISSDTD2020.roi_name(idx))
    end
%     pause
end
save(fullfile(s_dir, 'LI_glasser_bilateral_net_12.mat'), 'glass_net_regions', 'glass_net_L_idx','glass_net_R_idx', 'glass_net_R_label', 'glass_net_L_label');

%% Lateral ROIs for LI analysis (ROIs were provided by Joe)
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["parcel_num", "parcel_name", "num_vertices", "roi_num", "roi_name", "roi_count", "roi_percentage"];
opts.VariableTypes = ["double", "string", "double", "double", "string", "double", "double"];
opts = setvaropts(opts, 2, "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 5], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
glasserlateral = readtable("/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/glasser atlas/glasser_lateral.tsv", opts);
clear opts

glass_roi_L_idx = []; 
glass_roi_R_idx = [];
glass_roi_R_name = [];
glass_roi_L_name = [];
for i=1:length(glasserlateral.roi_name)
    if contains(glasserlateral.roi_name(i),'L_')
        glass_roi_L_idx = [glass_roi_L_idx, i]; 
        glass_roi_L_name = [glass_roi_L_name; glasserlateral.parcel_name(i)];
    elseif contains(glasserlateral.roi_name(i),'R_')
        glass_roi_R_idx = [glass_roi_R_idx, i];
        glass_roi_R_name = [glass_roi_R_name; glasserlateral.parcel_name(i)];
    end
end

s_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/glasser atlas';
glass_roi_name = [glass_roi_L_name; glass_roi_R_name];

%- Lookup
glass_roi_L_idx = []; 
glass_roi_R_idx = [];
for i = 1:length(glass_roi_L_name)
    if ~isempty(find(contains(glasserlateral.parcel_name, glass_roi_L_name(i))==1, 1))
        idx = find(contains(glasserlateral.parcel_name, glass_roi_L_name(i))==1);
        glass_roi_L_idx = [glass_roi_L_idx, idx];
        disp([glasserlateral.parcel_name(idx), glass_roi_L_name(i)])
    end
end
for i = 1:length(glass_roi_R_name)
    if ~isempty(find(contains(glasserlateral.parcel_name, glass_roi_R_name(i))==1, 1))
        idx = find(contains(glasserlateral.parcel_name, glass_roi_R_name(i))==1);
        glass_roi_R_idx = [glass_roi_R_idx, idx];
        disp([glasserlateral.parcel_name(idx), glass_roi_R_name(i)])
    end
end

save(fullfile(s_dir, 'LI_glasser_rois_48.mat'), 'glass_roi_name','glass_roi_R_idx','glass_roi_L_idx');

%% Read networks for LI analysis (Networks were provided by Joe)
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["parcel_num", "parcel_name", "num_vertices", "roi_num", "roi_name", "roi_count", "roi_percentage"];
opts.VariableTypes = ["double", "string", "double", "double", "string", "double", "double"];
opts = setvaropts(opts, 2, "WhitespaceRule", "preserve");
opts = setvaropts(opts, [2, 5], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
glasserROISSDTD2020 = readtable("/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/glasser atlas/glasser_ROIS_SDTD_2020.tsv", opts);
clear opts

glass_net_regions = unique(glasserROISSDTD2020.roi_name);


glass_net_L_idx = [];  glass_net_R_idx = [];
glass_net_L_label = [];  glass_net_R_label = [];
k1 = 1; k2 = 1;
for i=1:length(glass_net_regions)
    disp(glass_net_regions(i))
    if contains(glass_net_regions(i),'L')
%         glass_roi_L_idx = [glass_roi_L_idx, i];
        idx = find(contains(glasserROISSDTD2020.roi_name, glass_net_regions(i))==1);
        glass_net_L_idx{k1} = idx; 
        glass_net_L_label{k1} = glasserROISSDTD2020.parcel_name(idx);
        k1 = k1 + 1;
        disp(glasserROISSDTD2020.roi_name(idx))
    elseif contains(glass_net_regions(i),'R')
%         glass_roi_R_idx = [glass_roi_L_idx, i];
        idx = find(contains(glasserROISSDTD2020.roi_name, glass_net_regions(i))==1);
        glass_net_R_idx{k2} = idx;
        glass_net_R_label{k2} = glasserROISSDTD2020.parcel_name(idx);
        k2 = k2 + 1;
        disp(glasserROISSDTD2020.roi_name(idx))
    end
    pause
end
save(fullfile(s_dir, 'LI_glasser_net_12.mat'), 'glass_net_regions', 'glass_net_L_idx','glass_net_R_idx', 'glass_net_R_label', 'glass_net_L_label');

%% Read mask (fMRI_1)
cd('/group/jbinder/ECP/MEG/laterality_index')
ROI_AllCortex_Mask_new = ft_read_mri('lateral_roi.nii.gz');

cfg = [];
cfg.method        = 'ortho';
ft_sourceplot(cfg, ROI_AllCortex_Mask_new);

cfg = [];
cfg.funparameter = 'anatomy';
cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projmethod     = 'nearest'; 
cfg1 = ft_sourceplot(cfg, ROI_AllCortex_Mask_new);
view([-90,0])

%% Read mask (fMRI)
cd('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/ECP_masks')
ROI_AllCortex_Mask_new = ft_read_mri('lateral_roi.nii.gz');
MNI_ROIs_SDTD_2020 = ft_read_mri('MNI_ROIs_SDTD_2020.nii.gz');

cfg = [];
cfg.method        = 'ortho';
ft_sourceplot(cfg, ROI_AllCortex_Mask_new);
ft_sourceplot(cfg, MNI_ROIs_SDTD_2020);

cfg = [];
cfg.funparameter = 'anatomy';
cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projmethod     = 'nearest'; 
ft_sourceplot(cfg, MNI_ROIs_SDTD_2020);
view([-90,0])

cfg = [];
cfg.funparameter = 'anatomy';
cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projmethod     = 'nearest'; 
ft_sourceplot(cfg, ROI_AllCortex_Mask_new);
view([-90,0])

%% Read mask (fMRI)
cd('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/ECP_masks')

ROI_AllCortex_Mask = ft_read_mri('ROI_AllCortex_Mask.nii');
ROI_Angular_Mask = ft_read_mri('ROI_Angular_Mask.nii');
ROI_Temporal_Mask = ft_read_mri('ROI_Temporal_Mask.nii');
ROIs_NoOccipital = ft_read_mri('ROIs_NoOccipital.nii');
ROI_Frontal_Mask = ft_read_mri('ROI_Frontal_Mask.nii');

HalfBrain_ROI_Angular_Mask = ft_read_mri('HalfBrain_ROI_Angular_Mask.nii');
% atlas = ft_read_atlas('ROI_AllCortex_Mask.nii');

%% HCP nii
nii_hcp = '/data/MEG/Vahab/Github/MCW_MEGlab/FT/Atlas/HCP/HCP atlas for Brainstorm/MMP_in_MNI_symmetrical.nii';
hcp_nii = ft_read_mri(nii_hcp);

cfg = [];
cfg.method        = 'ortho';
ft_sourceplot(cfg, hcp_nii);

atlas_hcp = ft_read_atlas(nii_hcp);

cfg = [];
cfg.funparameter = 'anatomy';
cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projmethod     = 'nearest'; 
ft_sourceplot(cfg, hcp_nii);
view([-90,0])

cfg = []; cfg.resolution = 1; cfg.dim = [256 256 256];
hcp_nii1 = ft_volumereslice(cfg, hcp_nii); ft_sourceplot([], hcp_nii);

%%
close all
mask_in = ROI_Angular_Mask;
mask_in = ROI_Temporal_Mask;
mask_in = ROI_AllCortex_Mask;
mask_in = ROIs_NoOccipital;
mask_in = ROI_Frontal_Mask;

%%
cfg = []; cfg.resolution = 1; cfg.dim = [256 256 256];
mask_in = ft_volumereslice(cfg, mask_in); ft_sourceplot([], mask_in);

%%
mask = 'anatomy';
mask_in = mask_in; idx = mask_in.(mask) > 0; mask_in.(mask)(idx) = 1;
% idx = find(mask_in.(mask) > 0); idx_uni = unique(mask_in.(mask)(idx));

hcp_nii_mask = hcp_nii1; hcp_nii_mask.(mask) = hcp_nii_mask.(mask).*mask_in.(mask);
ft_sourceplot([], hcp_nii_mask);

% close all
% cfg = [];
% cfg.funparameter = 'anatomy';
% cfg.method = 'surface';
% cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% % cfg.surfinflated   = 'surface_inflated_left_caret.mat';
% cfg.projmethod     = 'nearest'; 
% ft_sourceplot(cfg, hcp_nii_mask);
% view([-90,0])

cfg = [];
cfg.mask = mask;
cfg.colormap = [];
cfg.thre = 0;
cfg.surfinflated = 'surface_inflated_both_caret.mat';
cfg.views = [-90,0];
cfg.tit = '';
cfg.saveflag = [];
do_mapvis(cfg, hcp_nii_mask)

idx = find(hcp_nii_mask.(mask) > 0); idx_uni = unique(hcp_nii_mask.(mask)(idx)); disp(idx_uni);

cfg = [];
cfg.atlas = atlas;
cfg.src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data_info/cortex_pial_low.fs';
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.lat_index = [idx_L, idx_R];
cfg.rois = rois;
do_plot_HCP_atlas(cfg);

%% AAL atlas
% cfg = [];
% cfg.funparameter = 'parcellation';
% cfg.method = 'surface';
% cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% cfg.projmethod     = 'nearest'; 
% ft_sourceplot(cfg, atlas);
% view([0,90])

cfg = [];
cfg.method        = 'ortho';
% ft_sourceplot(cfg, ROI_AllCortex_Mask);
% ft_sourceplot(cfg, ROI_Angular_Mask);
% ft_sourceplot(cfg, ROI_Temporal_Mask);
% ft_sourceplot(cfg, ROIs_NoOccipital);
ft_sourceplot(cfg, HalfBrain_ROI_Angular_Mask);

%%
mask = 'anatomy';
mask_val = (ROI_AllCortex_Mask.(mask) + ROI_Angular_Mask.(mask) ...
    + ROI_Temporal_Mask.(mask) + ...
    ROI_Frontal_Mask.(mask) + ROIs_NoOccipital.(mask));
% ./5;

mask_all = ROI_AllCortex_Mask;
mask_all.(mask) = mask_val;


cfg = [];
cfg.method        = 'ortho';
ft_sourceplot(cfg, mask_all);


cfg = [];
cfg.funparameter = 'anatomy';
cfg.method = 'surface';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% cfg.projmethod     = 'nearest'; 
ft_sourceplot(cfg, mask_all);
view([-90,0])



