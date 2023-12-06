function [LI_ROIval,pow] = do_LI_magnitude(cfg_main)
% Script for computing laterality indices from MEG language activation dSPM maps

% Load dSPM image grid and scout information
ImageGridAmp = abs(cfg_main.d_in);
sScout = cfg_main.atlas;

% Get left and right subregions from scout data
LHscout = [];
for i = 1:length(cfg_main.idx_L)
    LHscout = [LHscout, sScout.Scouts(cfg_main.idx_L(i)).Vertices];
end

RHscout = [];
for i = 1:length(cfg_main.idx_R)
    RHscout = [RHscout, sScout.Scouts(cfg_main.idx_R(i)).Vertices];
end

% Extract amplitude values for left and right subregions
LHvals = ImageGridAmp(LHscout);
RHvals = ImageGridAmp(RHscout);

% Calculate maximum values for left and right subregions
LH_max = max(LHvals(:));
RH_max = max(RHvals(:));
ROIMax = max(LH_max, RH_max);

% Set the threshold based on the chosen type
switch cfg_main.Threshtype
    case 1
        threshold = cfg_main.thre * max(ImageGridAmp(:)); % Global max threshold
    case 2
        threshold = cfg_main.thre * max(ImageGridAmp(:)); % Time max threshold
    case 3
        threshold = cfg_main.thre * ROIMax; % ROI max threshold
end

% Mag significant voxels in each hemisphere
pow_left = sum(ImageGridAmp(LHvals > threshold));
pow_right = sum(ImageGridAmp(RHvals > threshold));

% Calculate laterality index and total significant voxels
LI_ROIval = 100 * ((pow_left - pow_right) / (pow_left + pow_right));

pow = struct('left', pow_left, 'right', pow_right);


% Display results
% fprintf('Total Significant Voxels: %d\n', ROIcount);
% fprintf('Laterality Index: %.2f%%\n', LI_ROIcount);
