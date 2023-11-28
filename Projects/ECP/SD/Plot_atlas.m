
clear; clc, close('all'); warning off,

Run_load_surface_template

src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat';
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';

clc
cfg = []; 
cfg.src_fname = src_fname;
cfg.glass_dir = glass_dir;
cfg.glass_atlas = glass_atlas;
Data_hcp_atlas = ecpfunc_hcp_atlas2(cfg);

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
glasserlateral = readtable("/group/jbinder/ECP/MEG/laterality_index/bilateral_glasser_lateral.tsv", opts);
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

glass_roi_name = [glass_roi_L_name; glass_roi_R_name];

%- Lookup (sanity check of ROIs)
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

%% Lateral ROIs (similar to with fMRI)
glass_roi_lat_R_name = [];
glass_roi_lat_L_name = [];
for i=1:length(glass_roi_L_idx)
    idx_lateral_L(i) = find(contains(rois, glasserlateral.parcel_name(glass_roi_L_idx(i)))==1);
    glass_roi_lat_L_name = [glass_roi_lat_L_name; glasserlateral.parcel_name(glass_roi_L_idx(i))];
end
for i=1:length(glass_roi_R_idx)
    idx_lateral_R(i) = find(contains(rois, glasserlateral.parcel_name(glass_roi_R_idx(i)))==1);
    glass_roi_lat_R_name = [glass_roi_lat_R_name; glasserlateral.parcel_name(glass_roi_R_idx(i))];
end

% Plot HCP atlas
cfg = struct();
cfg.atlas = atlas;
cfg.src_fname = src_fname;
cfg.sel = 'roi'; % 'whole', 'left', 'right', 'roi';
cfg.index_L = idx_lateral_L';
cfg.index_R = idx_lateral_R';
cfg.rois = rois;
cfg.rois_sel = 1:length(glass_roi_L_idx);
cfg.title = '';
do_plot_HCP_atlas(cfg);

%%
s_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';
save(fullfile(s_dir, 'LI_glasser_lateral_rois.mat'), 'idx_lateral_L','idx_lateral_R', 'glass_roi_lat_L_name', 'glass_roi_lat_R_name');



