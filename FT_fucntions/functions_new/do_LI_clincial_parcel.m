function [LI_ROIcount, ROIcount] = do_LI_clincial_parcel(cfg_main)
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
LHvals = ImageGridAmp(LHscout(:),:);
RHvals = ImageGridAmp(RHscout(:),:);

%% PCA analysis
% COEFF = pca(LHvals'); pcaLHvals = LHvals' * COEFF(:,1); %figure, plot(pcaLHvals)
% COEFF = pca(RHvals'); pcaRHvals = RHvals' * COEFF(:,1); %figure, plot(pcaRHvals)
% 
% LHvals = pcaLHvals; RHvals = pcaRHvals;
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
end

% Count the number of significant voxels in each hemisphere
L_ROIcount = sum(LHvals(:) > threshold);
R_ROIcount = sum(RHvals(:) > threshold);

% Calculate laterality index and total significant voxels
LI_ROIcount = 100 * ((L_ROIcount - R_ROIcount) / (L_ROIcount + R_ROIcount));
ROIcount = L_ROIcount + R_ROIcount;

% Display results
% fprintf('Total Significant Voxels: %d\n', ROIcount);
% fprintf('Laterality Index: %.2f%%\n', LI_ROIcount);
