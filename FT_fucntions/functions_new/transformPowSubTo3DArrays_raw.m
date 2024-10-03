function LI_hc_pow_sub = transformPowSubTo3DArrays_raw(pow_sub)
    % Initialize the output structure
    LI_hc_pow_sub = struct('left_raw', [], 'right_raw', []);
    
    % Determine the sizes
    numSubjects = size(pow_sub, 1); % Number of subjects
    numIntervals = size(pow_sub, 2); % Number of intervals
    if ~isempty(pow_sub)
        numElements = size(pow_sub(1, 1).left, 1); % Length of the vector in each left/right field
    else
        numElements = 0; % Default in case pow_sub is empty
    end
    
    % Initialize the 3D arrays within the structure
    LI_hc_pow_sub.left = zeros(numSubjects, numIntervals, numElements);
    LI_hc_pow_sub.right = zeros(numSubjects, numIntervals, numElements);
    
    % Populate the 3D arrays from the struct array
    for i = 1:numSubjects
        for j = 1:numIntervals
            LI_hc_pow_sub.left(i, j, :) = pow_sub(i, j).left_raw;
            LI_hc_pow_sub.right(i, j, :) = pow_sub(i, j).right_raw;
        end
    end
end
