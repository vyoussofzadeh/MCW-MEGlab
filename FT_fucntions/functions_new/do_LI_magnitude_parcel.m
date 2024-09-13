function [LI_ROIval,pow] = do_LI_magnitude_parcel(cfg_main)
% Script for computing laterality indices from MEG language activation dSPM maps

% Load dSPM image grid and scout information
ImageGridAmp = abs(cfg_main.d_in);
% ImageGridAmp = (cfg_main.d_in);
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
LHvals = ImageGridAmp(LHscout,:);
RHvals = ImageGridAmp(RHscout,:);

%% PCA analysis
COEFF = pca(LHvals'); pcaLHvals = LHvals' * COEFF(:,1); %figure, plot(pcaLHvals)
COEFF = pca(RHvals'); pcaRHvals = RHvals' * COEFF(:,1); %figure, plot(pcaRHvals)

LHvals = pcaLHvals; RHvals = pcaRHvals;

%%

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
% pow_left = sum((LHvals(:) > threshold));
% pow_right = sum((RHvals(:) > threshold));

%%
% minPowerThreshold = 1;
minPowerThreshold = 0.1 * max(ImageGridAmp(:));  % Example threshold definition

% Apply minimum power threshold along with the existing maximum threshold
LHvals_filtered = LHvals(LHvals > threshold & LHvals > minPowerThreshold);
RHvals_filtered = RHvals(RHvals > threshold & RHvals > minPowerThreshold);

% % Now, calculate power using the filtered values
pow_left = sum(LHvals_filtered);
pow_right = sum(RHvals_filtered);

%%
% Using mean instead of sum for significant voxels in each hemisphere (DID NOT WORK!)
% pow_left = mean(ImageGridAmp(LHvals(:) > threshold));
% pow_right = mean(ImageGridAmp(RHvals(:) > threshold));

% normalize to the size of parcel
pow_left = pow_left/length(LHvals);
pow_right = pow_right/length(RHvals);


% Calculate laterality index and total significant voxels
LI_ROIval = 100 * ((pow_left - pow_right) / (pow_left + pow_right));

pow = struct('left', pow_left, 'right', pow_right);


% Display results
% fprintf('Total Significant Voxels: %d\n', ROIcount);
% fprintf('Laterality Index: %.2f%%\n', LI_ROIcount);


%% Assuming you define a minimum power threshold:
% minPowerThreshold = 0.4 * max(ImageGridAmp(:));  % Example threshold definition
% minPowerThreshold = 0;
% 
% % Apply the threshold filtering
% LHvals_filtered = LHvals(LHvals > threshold & LHvals > minPowerThreshold);
% RHvals_filtered = RHvals(RHvals > threshold & RHvals > minPowerThreshold);
% 
% % Calculate power using the filtered values
% pow_left_filtered = sum(LHvals_filtered);
% pow_right_filtered = sum(RHvals_filtered);
% 
% % Normalize to the size of parcels, ensuring only counts that meet both thresholds are considered
% normalized_pow_left = pow_left_filtered / length(LHvals_filtered);
% normalized_pow_right = pow_right_filtered / length(RHvals_filtered);
% 
% % Calculate the Laterality Index using the normalized powers
% LI_ROIval = 100 * ((normalized_pow_left - normalized_pow_right) / (normalized_pow_left + normalized_pow_right));
% 
% pow = struct('left', normalized_pow_left, 'right', normalized_pow_right);


%%
% % Assuming you define a minimum power threshold:
% minPowerThreshold = cfg_main.minPowerThreshold;  % Define minimum threshold as part of your configuration
% 
% % Extract amplitude values for left and right subregions that exceed the minimum threshold
% LHvals = ImageGridAmp(LHscout,:);
% RHvals = ImageGridAmp(RHscout,:);
% 
% % Apply minimum power threshold along with the existing maximum threshold
% LHvals_filtered = LHvals(LHvals > threshold & LHvals > minPowerThreshold);
% RHvals_filtered = RHvals(RHvals > threshold & RHvals > minPowerThreshold);
% 
% % Now, calculate power using the filtered values
% pow_left = sum(LHvals_filtered);
% pow_right = sum(RHvals_filtered);
% 
% % Normalize to the size of parcels, but ensure you're using counts that meet both thresholds
% normalized_pow_left = pow_left / numel(LHvals_filtered);
% normalized_pow_right = pow_right / numel(RHvals_filtered);
% 
% % Calculate laterality index using normalized powers
% LI_ROIval = 100 * ((normalized_pow_left - normalized_pow_right) / (normalized_pow_left + normalized_pow_right));
% 
% pow = struct('left', normalized_pow_left, 'right', normalized_pow_right);

