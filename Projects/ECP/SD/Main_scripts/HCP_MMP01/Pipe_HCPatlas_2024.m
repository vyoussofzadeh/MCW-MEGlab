clear; clc, close('all'); warning off,

%% Paths
% restoredefaultpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/run')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/function')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/tools/helpful_tools/daviolinplot/daboxplot')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Main_scripts/run')
Run_setpath

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other/')

src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim';
% glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat';
glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Main_scripts/HCP_MMP01/scout_HCP_MMP1_360.mat';

glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';

LI_glasser_lateral_rois = load(fullfile(glass_dir,'LI_glasser_lateral_rois.mat'));

BS_atlas = fullfile(glass_dir,'LI_glasser_lateral_rois.mat');

cfg = struct('src_fname', src_fname, 'glass_dir', glass_dir, 'glass_atlas', glass_atlas, 'plotflag', 0, 'BS_atlas', BS_atlas);
Data_hcp_atlas = ecpfunc_hcp_atlas4(cfg);


save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Main_scripts/HCP_MMP01';

cfg = []; cfg.glass_dir = save_dir; 
cfg.Data_hcp_atlas = Data_hcp_atlas; 
ecpfunc_hcp_atlas_save(cfg);

%%
close all
cfg = [];
cfg.src_fname = src_fname;
cfg.network_sel = [1,2,6];
cfg.network_sel = [4];
cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.plotflag = 1;
cfg.fixedcolor = [0,0.7,0];
[idx_L, idx_R, src]  = do_plot_hcp_network2(cfg);
net_rois = 'ftp';

disp(Data_hcp_atlas.groups_labels)


cfg = []; cfg.idx_L = idx_L; cfg.idx_R = idx_R; cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.src = src;

save_dir_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/atlas_roi';
outdir = save_dir_atlas;

cfg.export = 0; cfg.savedir = fullfile(outdir,'group');
cfg.network_sel = [1,2,6]; 
cfg.plotflag = 0;
do_plot_hcp_network2(cfg);
net_label = 'Fronto_tempro_pri';
cfg.network_sel = [11]; do_map_HCP_net_sel2(cfg);

colorcode = {[0,114,189]; [217,83,25]; [237,177,32];[126,47,142]; [119,172,48]};
network_sel = [1,2,5,6,11]; % LI networks compared between MEG and fMRI

% Inspecting atlas networks common between MEG and fMRI
for i = 1:length(colorcode)
    cfg.network_sel = network_sel(i);
    cfg.fixedcolor = colorcode{i}/256;
    do_map_HCP_net_sel2(cfg);
    title(Data_hcp_atlas.groups_labels{cfg.network_sel});
    
    cfg2 = [];
    cfg2.outdir = save_dir_atlas;
    filename = Data_hcp_atlas.groups_labels{cfg.network_sel};
    cfg2.filename = filename;
    cfg2.type = 'fig';
    %         do_export_fig(cfg2)
end

%% Plotting using Conn toolbox
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/brewermap')

% SPM
spmpath = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/SPM/spm12_2021/spm12';
addpath(spmpath)
spm_get_defaults

% Conn
conn_path = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Conn/conn';
addpath(genpath(conn_path));

% look up rois

% Open the text file
fid = fopen('MMP_in_MNI_symmetrical.txt', 'r');

% Read each line as a string
rois_data = textscan(fid, '%s', 'Delimiter', '\n');

% Close the file
fclose(fid);

% Extract the cell array of ROIs
rois = rois_data{1};

% Assume rois is your 180x1 cell array described above.

n = numel(rois);
roiIndex = zeros(n, 1);
roiNames = cell(n, 1);

for i = 1:n
    thisLine = rois{i};              % e.g. '1 L_V1_ROI'
    tokens = strsplit(strtrim(thisLine));
    roiIndex(i) = str2double(tokens{1});  % convert '1' -> 1
    roiNames{i} = tokens{2};             % 'L_V1_ROI'
end

rois_id = [1,2,6,11]; idx_rois = [];

for j=1:length(rois_id)
    
    labels = Data_hcp_atlas.glass_net_L_label{rois_id(j)};
    idx = [];
    for i = 1:length(labels)
        
        thislabel = labels{i};
        thislabel = regexprep(thislabel, '\sL$', '');

        % Compare against all labels
        idx(i) = find(strcmp(roiNames, thislabel));
        
        if ~isempty(idx)
            fprintf('ROI %s found at index %d\n', thislabel, idx(i));
        else
            fprintf('ROI %s NOT found in the atlas.\n', thislabel);
        end
    end
    idx_rois{j} = idx;
end

disp('ang:')
disp(idx_rois{1})

disp('front:')
disp(idx_rois{2})

disp('temp:')
disp(idx_rois{3})

disp('lateral:')
disp(idx_rois{4})

cd('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/neurovault_MMP/MMP 1.0 MNI projections')
conn_mesh_display('MMP_in_MNI_symmetrical.nii');


hcp_atlas_ft = ft_read_mri('MMP_in_MNI_symmetrical.nii');

unique(hcp_atlas_ft.anatomy)


%% getting color match with new HCP
hcp_old = load('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/scout_mmp_in_mni_corr_updated.nii_362.mat');
hcp_new = load('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/scout_HCP_MMP1_360.mat');

%
% We'll skip the last 2 Scouts if you want (per your example).
hcp_old_label = {hcp_old.Scouts(1:end-2).Label}; 
hcp_new_label = {hcp_new.Scouts.Label(1:end-3)};

hcp_old_label = hcp_old_label(:);  % convert to column cell array
hcp_new_label = hcp_new_label(:);  % convert to column cell array

% 1. Remove the trailing space+L from each label
hcp_new_label_fixed = regexprep(hcp_new_label, '\sL$', '');
hcp_new_label_fixed = regexprep(hcp_new_label_fixed, '\sR$', '');


[foundMask, idxMatch] = ismember(hcp_old_label, hcp_new_label_fixed);

% foundMask is a logical array of size size(hcp_old_label)
% idxMatch gives the matching index in hcp_new_label, or 0 if not found.

% For debugging, see which ones matched vs. not:
matchedOldLabels = hcp_old_label(foundMask);
matchedIdx       = idxMatch(foundMask);

notFoundLabels   = hcp_old_label(~foundMask);

disp('Matched:');
disp([matchedOldLabels(:), num2cell(matchedIdx(:))]);

disp('Not found:');
disp(notFoundLabels);

%% 1) Load the old and new HCP data
cd('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm')
hcp_old = load('scout_mmp_in_mni_corr_updated.nii_362.mat');
hcp_new = load('scout_HCP_MMP1_360.mat');

%% 2) Extract label arrays
% Suppose hcp_old has >=180 Scouts, hcp_new has 360
oldLabels = {hcp_old.Scouts.Label};    % e.g. Nx1 or 1xN cell array
newLabels = {hcp_new.Scouts.Label};    % e.g. 360x1 or 1x360

% If newLabels has " L" at the end, remove it
newLabelsFixed = regexprep(newLabels, '\sL$', '');  % "L_V1_ROI L" -> "L_V1_ROI"
newLabelsFixed = regexprep(newLabelsFixed, '\sR$', '');

%% 3) Update the color for the first 180 scouts in hcp_old
numScoutsToUpdate = 360;
colorMapMatrix = zeros(numScoutsToUpdate, 3);  % to store the 180 color values

for i = 1:numScoutsToUpdate
    % The label in hcp_old
    thisOldLabel = oldLabels{i};   
    
    % Find matching label index in hcp_new
    idxMatch = find(strcmp(newLabelsFixed, thisOldLabel));
    if ~isempty(idxMatch)
        % If a match is found, update hcp_old color from hcp_new
        hcp_old.Scouts(i).Color = hcp_new.Scouts(idxMatch).Color;
        
        % Store the color in our 180x3 matrix
        colorMapMatrix(i, :) = hcp_new.Scouts(idxMatch).Color;
    else
        % No match found, optionally print a warning or skip
        fprintf('Warning: No match found for label "%s".\n', thisOldLabel);
    end
end

%% 4) Check the colorMapMatrix for the updated colors
% For example, look at the first 5
hcp_color = colorMapMatrix(1:180,:);

%% 5) Save the updated hcp_old (and the colorMapMatrix) to a new .mat file
save('hcp_color_for_conn.mat', 'hcp_color');
save('scout_mmp_in_mni_corr_updated.nii_362_newColor.mat','-struct', 'hcp_old')





