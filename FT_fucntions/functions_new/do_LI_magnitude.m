function [LI_ROIval,pow] = do_LI_magnitude(cfg_main)
% Script for computing laterality indices from MEG language activation dSPM maps

% Load dSPM image grid and scout information
ImageGridAmp = abs(cfg_main.d_in);
% ImageGridAmp = 10.^(ImageGridAmp / 10);

% ImageGridAmp = (cfg_main.d_in);
sScout = cfg_main.atlas;

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
        %         minPowerThreshold =  0.01 * cfg_main.globalmax; % ROI max threshold
end

%%
% Mag significant voxels in each hemisphere
pow_left = sum(ImageGridAmp(LHvals(:) > threshold));
pow_right = sum(ImageGridAmp(RHvals(:) > threshold));

pow_left = pow_left/length(LHvals);
pow_right = pow_right/length(RHvals);

% pow_left = median(ImageGridAmp(LHvals(:) > threshold));
% pow_right = median(ImageGridAmp(RHvals(:) > threshold));

%%
% % minPowerThreshold = 1;
% minPowerThreshold = 0.1 * max(ImageGridAmp(:));  % Example threshold definition
%
% % Apply minimum power threshold along with the existing maximum threshold
% LHvals_filtered = LHvals(LHvals > threshold & LHvals > minPowerThreshold);
% RHvals_filtered = RHvals(RHvals > threshold & RHvals > minPowerThreshold);
%
% % % Now, calculate power using the filtered values
% pow_left = sum(LHvals_filtered);
% pow_right = sum(RHvals_filtered);
%
% % normalize to the size of parcel
% pow_left = pow_left/length(LHvals_filtered);
% pow_right = pow_right/length(RHvals_filtered);

%%
% pow_left = mean(LHvals_filtered);
% pow_right = mean(RHvals_filtered);

%%
% Calculate laterality index and total significant voxels
LI_ROIval = 100 * ((pow_left - pow_right) / (pow_left + pow_right));

pow = struct('left', pow_left, 'right', pow_right);


% Display results
% fprintf('Total Significant Voxels: %d\n', ROIcount);
% fprintf('Laterality Index: %.2f%%\n', LI_ROIcount);
