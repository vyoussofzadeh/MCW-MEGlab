% Basic consistency checks
assert(size(VolResults.GridLoc, 1) == ...
       size(VolResults.ImageGridAmp, 1), ...
       'GridLoc and ImageGridAmp have inconsistent dimensions.');

assert(size(CortexMat.Vertices, 2) == 3, ...
       'Cortex vertices must be an N-by-3 matrix.');

% Shepard inverse-distance interpolation:
% destination = cortical vertices
% source      = volume grid coordinates
nNeighbors = 8;
excludeDistance = -10;  % Ignore volume points >10 mm from cortex
distanceExponent = 2;

W = bst_shepards( ...
    CortexMat.Vertices, ...
    VolResults.GridLoc, ...
    nNeighbors, ...
    excludeDistance, ...
    distanceExponent);

fprintf('Interpolation matrix: %d cortex vertices x %d volume points\n', ...
    size(W,1), size(W,2));

% Interpolate all available time points
SurfaceValues = W * VolResults.ImageGridAmp;

%%
SurfResults = SurfTemplate;

% Replace template activity with interpolated DICS values
SurfResults.ImageGridAmp = SurfaceValues;
SurfResults.Time         = VolResults.Time;
SurfResults.Comment      = [VolResults.Comment ' | volume-to-surface'];
SurfResults.Function     = VolResults.Function;
SurfResults.DataFile     = VolResults.DataFile;
SurfResults.nComponents  = 1;

% Surface-source files do not use explicit GridLoc
SurfResults.GridLoc = [];

if isfield(SurfResults, 'GridOrient')
    SurfResults.GridOrient = [];
end

% Preserve relevant provenance
SurfResults = bst_history('add', SurfResults, 'compute', ...
    'Interpolated from volume grid to cortex with bst_shepards (8 neighbors, 10-mm limit)');