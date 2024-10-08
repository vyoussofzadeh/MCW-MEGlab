function [LI_ROIval,pow] = do_LI_magnitude(cfg_main)
% Script for computing laterality indices from MEG language activation dSPM maps

% Load dSPM image grid and scout information
ImageGridAmp = cfg_main.d_in;

sScout = cfg_main.atlas;

if size(ImageGridAmp,1) > 360
    cfg_main.parcellaion = 0;
end

if cfg_main.parcellaion == 1
    
    % Extract amplitude values for left and right subregions
    LHvals = ImageGridAmp(cfg_main.idx_L,:);
    RHvals = ImageGridAmp(cfg_main.idx_R,:);
else
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
    LHvals = ImageGridAmp(LHscout,:);
    RHvals = ImageGridAmp(RHscout,:);
end

if cfg_main.applymean == 1
    % Mag significant voxels in each hemisphere
    LHvals = LHvals(1,:);
    RHvals = RHvals(1,:);
end

% Calculate maximum values for left and right subregions
LH_max = max(LHvals(:));
RH_max = max(RHvals(:));
ROIMax = max(LH_max, RH_max);

% Set the threshold based on the chosen type
switch cfg_main.Threshtype
    case 1
        threshold = cfg_main.thre * cfg_main.globalmax; % Global max threshold
    case 2
        threshold = cfg_main.thre * max(ImageGridAmp(:)); % Time max threshold
    case 3
        threshold = cfg_main.thre * ROIMax; % ROI max threshold
    case 4
        threshold = cfg_main.thre * cfg_main.globalwindowthresh;
end

%% raw time-series, aka raw-SNR
pow_left_raw  = mean(LHvals(:)); 
pow_right_raw  = mean(RHvals(:));

%% thresholded timeseries, aka, response-SNR
if cfg_main.applymean == 1

    idx_left = LHvals > threshold; 
    pow_left = sum(LHvals(idx_left));
    
    idx_right = RHvals > threshold; 
    pow_right = sum(RHvals(idx_right));
    
else
    pow_left  = sum(LHvals(LHvals(:) > threshold)); 
    pow_left = pow_left/size(LHvals,1);
    
    pow_right = sum(RHvals(RHvals(:) > threshold)); 
    pow_right = pow_right/size(RHvals,1);

end

% Calculate laterality index and total significant voxels
LI_ROIval = 100 * ((pow_left - pow_right) / (pow_left + pow_right));

pow = struct('left', pow_left, 'right', pow_right, 'left_raw', pow_left_raw, 'right_raw', pow_right_raw);


% Display results
% fprintf('Total Significant Voxels: %d\n', ROIcount);
% fprintf('Laterality Index: %.2f%%\n', LI_ROIcount);
