function ecpfunc_hcp_atlas_save(cfg_main)

Data_hcp_atlas = cfg_main.Data_hcp_atlas;
glass_dir = cfg_main.glass_dir;

% Assuming Data_hcp_atlas is already computed and available
atlas = cfg_main.Data_hcp_atlas.atlas;

% Pre-allocate cell array for all new scouts
new_scouts_all = cell(1, length(Data_hcp_atlas.groups_labels));

% Iterate over each group
for group_idx = 1:length(Data_hcp_atlas.groups_labels)
    
    % Concatenate left and right subregions
    left_subregions = Data_hcp_atlas.glass_net_L_label{group_idx};
    right_subregions = Data_hcp_atlas.glass_net_R_label{group_idx};
    all_subregions = [left_subregions; right_subregions];
    
    % Get the list of all scout labels in the atlas
    rois = {atlas.Scouts.Label};
    
    % Find indices of matching scouts
    idx = zeros(1, length(all_subregions)); % Pre-allocate for performance
    for j = 1:length(all_subregions)
        match = find(strcmp(rois, all_subregions{j}), 1);
        if ~isempty(match)
            idx(j) = match;
        else
            idx(j) = NaN; % Assign NaN for no match to avoid indexing errors later
        end
    end
    
    % Filter out any NaN values from idx to avoid indexing errors
    valid_indices = idx(~isnan(idx));
    
    % Store only the valid scouts in the array
    if ~isempty(valid_indices)
        new_scouts_all{group_idx} = atlas.Scouts(valid_indices);
    else
        new_scouts_all{group_idx} = [];
    end
end

% Save each group's modified atlas to a new file
for i = 1:length(Data_hcp_atlas.groups_labels)
    
    if isempty(new_scouts_all{i})
        continue; % Skip saving if no scouts were found
    end
    
    % Prepare the filename
    new_atlas_filename = fullfile(glass_dir, ['Scout_', Data_hcp_atlas.groups_labels{i}, '.mat']);
    
    % Update the atlas with new scouts for the current group
    atlas.Scouts = new_scouts_all{i}; % Update the Scouts field
    atlas.Name = Data_hcp_atlas.groups_labels{i};
    
    % Save the modified atlas
    save(new_atlas_filename, '-struct','atlas'); % Removed '-struct' for correct saving of the full struct
    
    % Display the saving message
    disp(['Scout_atlas saved to: ', new_atlas_filename]);
    
end

cd(glass_dir)