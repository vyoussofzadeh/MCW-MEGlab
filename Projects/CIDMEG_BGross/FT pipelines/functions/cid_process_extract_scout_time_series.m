function scout_time_series = cid_process_extract_scout_time_series(sResults, sScout, method)
    % Extract scout-based time series from source results
    %
    % Input:
    %   sResults - Structure from in_bst_results containing source data
    %   sScout   - Brainnetome scout structure with region definitions
    %   method   - String specifying extraction method ('mean' or 'all')
    %
    % Output:
    %   scout_time_series - A structure with scout labels as fields and
    %                       corresponding time series as values.

    % Default method to 'mean' if not provided
    if nargin < 3
        method = 'mean';
    end
    
    % Check if the source data is available
    if ~isfield(sResults, 'ImageGridAmp') || isempty(sResults.ImageGridAmp)
        error('Source results do not contain ImageGridAmp field.');
    end
    
    % Initialize output structure
    scout_time_series = struct();
    
    % Loop through each scout in the atlas
    for iScout = 1:length(sScout.Scouts)
        % Get scout name and sanitize for use as a field name
        original_scout_name = sScout.Scouts(iScout).Label;
        sanitized_scout_name = matlab.lang.makeValidName(original_scout_name);
        
        % Get scout vertices
        scout_indices = sScout.Scouts(iScout).Vertices;
        
        % Ensure the vertices exist in the source grid
        if any(scout_indices > size(sResults.ImageGridAmp, 1))
            warning(['Vertices for scout "', original_scout_name, '" exceed the source grid size. Skipping.']);
            continue;
        end
        
        % Extract time series for the scout based on the chosen method
        switch method
            case 'mean'
                % Average across vertices
                time_series = mean(sResults.ImageGridAmp(scout_indices, :), 1);
            case 'all'
                % Keep all vertex time series
                time_series = sResults.ImageGridAmp(scout_indices, :);
            otherwise
                error('Invalid method. Choose "mean" or "all".');
        end
        
        % Save time series into the structure
        scout_time_series.(sanitized_scout_name) = time_series;
    end
    
    % Output log
    disp(['Extracted time series for ', num2str(length(fields(scout_time_series))), ' scouts.']);
end

