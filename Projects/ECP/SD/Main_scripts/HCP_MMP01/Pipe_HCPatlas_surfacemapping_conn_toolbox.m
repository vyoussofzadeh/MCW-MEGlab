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

%% Plotting using Conn toolbox
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/brewermap')
cd('/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/neurovault_MMP/MMP 1.0 MNI projections')

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

conn_mesh_display('MMP_in_MNI_symmetrical.nii');

hcp_atlas_ft = ft_read_mri('MMP_in_MNI_symmetrical.nii');
unique(hcp_atlas_ft.anatomy)

%%
ang = [16    17   141   143   145   150   151   152];
front = [72    65    88    91    92    32    58    74    75    84    76    66    94    12    68    67    70    63  73    86    87    69    71   111   169    79    80    82    81    26    89   179    77    85    62    97 170   171    83   165    98];
temp = [18   118   119   120   122   123   125   126   127   128   129   130   131   132   133   134   135   136 137   138   155   163   172   176   177];
lateral = [11    16    17    46    50    66    67    73    74    75    76    77    79    80    81    82    83    92 94    97   111   112   123   125   128   129   130   131   132   133   134   136   137   138   141   143 145   146   150   151   169   171   172   176   177];

ang_frontal_temp =  [16    17   141   143   145   150   151   152 ...
72    65    88    91    92    32    58    74    75    84    76    66   ...
94    12    68    67    70    63  73    86    87    69    71   111  ...
169    79    80    82    81    26    89   179    77    85   ...
62    97 170   171    83   165    98 18   118   119   120   122   123   ...
125   126   127   128   129   130   131   132   133   134   135   136 137   138   155   163   172   176   177];


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

%% --------------------------------------------------------------
% Define your four fixed colours (0255 ? 01)
colAngular  = [  0 114 188]/255;     % teal
colFrontal  = [216  83  25]/255;     % orange
colTemporal = [127  96   0]/255;     % brown
colLateral  = [126  47 141]/255;     % purple

% --------------------------------------------------------------
% ROI index lists (1180 only)
idxAngular  = idx_rois{1};
idxFrontal  = idx_rois{2};
idxTemporal = idx_rois{3};
idxLateral  = idx_rois{4};

% --------------------------------------------------------------
% 1) Colour-map for **lateral only**
cmapLat = zeros(180,3);                        % start black
cmapLat(idxLateral,:) = repmat(colLateral, numel(idxLateral), 1);

% --------------------------------------------------------------
% 2) Colour-map for **combined other regions**
cmapOther = zeros(180,3);
cmapOther(idxAngular ,:) = repmat(colAngular , numel(idxAngular ), 1);
cmapOther(idxFrontal ,:) = repmat(colFrontal , numel(idxFrontal ), 1);
cmapOther(idxTemporal,:) = repmat(colTemporal, numel(idxTemporal), 1);
% Lateral rows stay black ? easy to spot absence

% --------------------------------------------------------------
% (Optional) update atlas colours in-memory
paint = @(idx,rgb) arrayfun(@(ii)setfield(hcp_old.Scouts(ii),'Color',rgb),idx);

% if you want the atlas itself to show these same colours:
paint(idxAngular , colAngular );
paint(idxFrontal , colFrontal );
paint(idxTemporal, colTemporal);
paint(idxLateral , colLateral );

% --------------------------------------------------------------
% Save colour maps for CONN / custom scripts
atlasDir = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm';
cd(atlasDir)

save('hcp_cmap_lateral.mat' , 'cmapLat');
save('hcp_cmap_other.mat'   , 'cmapOther');

%% fixed colour
% [0,114,188] -> angular
% [216,83,25] -> frontal
% [126,47,141] -> lateral
% [127,96,0] -> temporal


%%
view([-90,0])
view([180,0])
view([180,-90])
view([90,0])
view([0,90])

view([-120,20])




