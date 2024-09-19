%% The CIDMEG project

% Export Association network ROIs into Brainstorm format (.mat)
% Written by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Latest update: 08/26/2024

%%
% Initial setup: Clear workspace, load scout data, and define regions
clear; clc;

% Define the directory for saving data
export_dir = '/group/bgross/work/CIDMEG/analysis/Pipelines/atlas_associationROIs';

% sScout = load('scout_Desikan-Killiany_68.mat');
sScout = load('scout_Desikan-Killiany_68_BS2024.mat');


% expectedRegions = {'superiorparietal', 'inferiorparietal', 'supramarginal', ...
%     'postcentral', 'precuneus', 'superiortemporal', 'middletemporal', ...
%     'inferiortemporal', 'bankssts', 'fusiform', 'transversetemporal', ...
%     'entorhinal', 'temporalpole', 'parahippocampal', 'lateraloccipital', ...
%     'lingual', 'cuneus', 'pericalcarine', 'posteriorcingulate', 'isthmuscingulate'};

expectedRegions = {'superiorparietal', 'inferiorparietal', 'supramarginal', ...
    'precuneus', 'superiortemporal', 'middletemporal', ...
    'inferiortemporal', 'bankssts', 'fusiform', 'transversetemporal', ...
    'entorhinal', 'temporalpole', 'parahippocampal', ...
    'posteriorcingulate', 'isthmuscingulate'};


% Normalize labels to ensure consistency and ease of comparison
normalizedRegions = cellfun(@(x) regexprep(x, '[ L|R]$', ''), {sScout.Scouts.Label}, 'UniformOutput', false);
normalizedRegions = cellfun(@(x) strtrim(x), normalizedRegions, 'UniformOutput', false);

% Find indices corresponding to the expected regions
selectedIndices = arrayfun(@(x) find(strcmp(normalizedRegions, expectedRegions{x})), 1:length(expectedRegions), 'UniformOutput', false);

% Display the indices for verification
disp('Indices of selected ROIs:');
for i = 1:length(expectedRegions)
    if isempty(selectedIndices{i})
        fprintf('%s: None found\n', expectedRegions{i});
    else
        fprintf('%s: %s\n', expectedRegions{i}, mat2str(selectedIndices{i}));
    end
end

% Convert selected indices from a cell array to a flat vector
allIndices = cell2mat(selectedIndices(:));
flatIndices = reshape(allIndices.', 1, []);  % Ensure a single row of indices

%% SAVE ROIS in BS format
% Verify and create the directory if it doesn't exist
if ~exist(export_dir, 'dir')
    fprintf('Creating directory: %s\n', export_dir);
    mkdir(export_dir);
end

% Update the scouts and save the new atlas
newScoutsAll = sScout.Scouts(flatIndices);
newatlas = sScout;
newatlas.Scouts = newScoutsAll;
newatlas.Name = 'association';
new_atlas_filename = fullfile(export_dir, 'scout_associationROIs_BS24.mat');

% Save the updated atlas structure
fprintf('Saving new atlas structure to %s\n', new_atlas_filename);
save(new_atlas_filename, '-struct', 'newatlas');
fprintf('Atlas saved successfully.\n');

