function LI_hc_pow_sub = transformCountSubTo3DArrays(count_sub)
    % Initialize the output structure
    LI_hc_pow_sub = struct('left', [], 'right', []);
    
    % Determine the sizes
    numSubjects = size(count_sub, 1); % Number of subjects
    numIntervals = size(count_sub, 2); % Number of intervals
    if ~isempty(count_sub)
        numElements = size(count_sub(1, 1).left, 1); % Length of the vector in each left/right field
    else
        numElements = 0; % Default in case count_sub is empty
    end
    
    % Initialize the 3D arrays within the structure
    LI_hc_count_sub.left = zeros(numSubjects, numIntervals, numElements);
    LI_hc_count_sub.right = zeros(numSubjects, numIntervals, numElements);
    
    % Populate the 3D arrays from the struct array
    for i = 1:numSubjects
        for j = 1:numIntervals
            LI_hc_count_sub.left(i, j, :) = count_sub(i, j).left;
            LI_hc_count_sub.right(i, j, :) = count_sub(i, j).right;
        end
    end
end
