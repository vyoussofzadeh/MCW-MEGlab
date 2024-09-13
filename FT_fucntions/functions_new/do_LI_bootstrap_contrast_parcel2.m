function [weighted_li, LI_max, num_threshvals] = do_LI_bootstrap_contrast_parcel2(cfg_main)

% MIN_NUM_THRESH_VOXELS is a predefined minimum number ...
% of voxels that are required for the computation to be considered valid.

wi = cfg_main.wi;
atlas = cfg_main.atlas;
idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;
divs = cfg_main.divs; % divs: number of divisions between 0 and max value in img
n_resampling = cfg_main.n_resampling; % 200
dvd = cfg_main.dvd; % 5
doavg = cfg_main.doavg;
downsamplerate = cfg_main.downsamplerate;

sinput = cfg_main.sinput;

RESAMPLE_RATIO = cfg_main.RESAMPLE_RATIO; % 0.75
MIN_NUM_THRESH_VOXELS = round(dvd / RESAMPLE_RATIO);

%%
% tmp = load(fullfile(cfg_main.BS_data_dir, sinput));

tmp_1 = load(fullfile(cfg_main.BS_data_dir, sinput{1}));
tmp_2 = load(fullfile(cfg_main.BS_data_dir, sinput{2}));

%%
sScout = atlas;

%% Parcellation (active)
% Get left and right subregions from scout data
LHscout = [];
for i = 1:length(idx_L)
    LHscout = [LHscout, sScout.Scouts(idx_L(i)).Vertices];
end

RHscout = [];
for i = 1:length(idx_R)
    RHscout = [RHscout, sScout.Scouts(idx_R(i)).Vertices];
end

%% Source timeseries and PCA analysis: Active condition
% Extract amplitude values for left and right subregions
LHvals = tmp_1.ImageGridAmp(LHscout(:),:);
RHvals = tmp_1.ImageGridAmp(RHscout(:),:);

% PCA analysis
COEFF = pca(LHvals'); pcaLHvals = LHvals' * COEFF(:,1); %figure, plot(pcaLHvals)
COEFF = pca(RHvals'); pcaRHvals = RHvals' * COEFF(:,1); %figure, plot(pcaRHvals)

LHval_act = pcaLHvals; RHval_act = pcaRHvals;

%% Source timeseries and PCA analysis: Baseline condition
% Extract amplitude values for left and right subregions
LHvals = tmp_2.ImageGridAmp(LHscout(:),:);
RHvals = tmp_2.ImageGridAmp(RHscout(:),:);

% PCA analysis
COEFF = pca(LHvals'); pcaLHvals = LHvals' * COEFF(:,1); %figure, plot(pcaLHvals)
COEFF = pca(RHvals'); pcaRHvals = RHvals' * COEFF(:,1); %figure, plot(pcaRHvals)

LHval_bsl = pcaLHvals; RHval_bsl = pcaRHvals;

% figure, plot(tmp_2.Time, LHval_bsl)
% figure, plot(tmp_1.Time, LHval_act)

%% Subtraction
LHval = LHval_act - LHval_bsl;
RLHval = RHval_act - RHval_bsl;

sTime = tmp_1.Time;

weighted_li = [];
% ft_progress('init', 'text',     'please wait ...');

for j = 1:size(wi, 1)
    
    %     ft_progress(j/size(wi, 1), 'Processing interval %d from %d', j, size(wi, 1));
    timind1 = nearest(sTime, wi(j,1));
    timind2 = nearest(sTime, wi(j,2));
    
    LHvals_interval = LHvals(timind1:timind2);
    RHvals_interval = RHvals(timind1:timind2);
    
    % downsample data
    LHvals_interval = downsample(LHvals_interval,downsamplerate);
    RHvals_interval = downsample(RHvals_interval,downsamplerate);

   % Calculate maximum values for left and right subregions
    LH_max = max(LHvals_interval);
    RH_max = max(RHvals_interval);
    ROIMax = max(LH_max, RH_max);
    
    %     divs = 5; % divs: number of divisions between 0 and max value in img
    lvals_nonnegative = LHvals_interval(LHvals_interval(:) >= 0);
    rvals_nonnegative = RHvals_interval(RHvals_interval(:) >= 0);
    
    % creating an array of values between 0 and the maximum amplitude
    % (ROIMax) observed in the brain regions of interest.
    % This array is linearly spaced and divided into divs number of divisions.
    % This creates a series of thresholds that will be used to progressively ...
    % include only those voxels with amplitude values above each threshold.
    threshvals = (0:(divs-1)) * (ROIMax / (divs - 1));
    
    lr_threshvals = cell(1, numel(threshvals));
    num_threshvals = 0;
    
    for i = 1:numel(threshvals)
        threshval = threshvals(i);
        l_threshvals = lvals_nonnegative(lvals_nonnegative(:) >= threshval);
        r_threshvals = rvals_nonnegative(rvals_nonnegative(:) >= threshval);
        
        %- Check if both hemispheres have enough voxels above the threshold
        % (MIN_NUM_THRESH_VOXELS)
        if (numel(l_threshvals) >= MIN_NUM_THRESH_VOXELS && numel(r_threshvals) >= MIN_NUM_THRESH_VOXELS)
            num_threshvals = num_threshvals + 1;
            lr_threshvals{num_threshvals} = {l_threshvals, r_threshvals};
        else
            break;
        end
    end
    interval_length(j) = length(lr_threshvals);
    
    TB_LIs = zeros(num_threshvals, n_resampling);
    
    for i = 1:num_threshvals
        l_set = lr_threshvals{i}{1};
        r_set = lr_threshvals{i}{2};
        
        l_n = max(round(numel(l_set) * RESAMPLE_RATIO), MIN_NUM_THRESH_VOXELS);
        r_n = max(round(numel(r_set) * RESAMPLE_RATIO), MIN_NUM_THRESH_VOXELS);
        
        l_indices = randi(numel(l_set), n_resampling, l_n);
        lactivity = sum(l_set(l_indices), 2);
        
        r_indices = randi(numel(r_set), n_resampling, r_n);
        ractivity = sum(r_set(r_indices), 2);
        
        TB_LIs(i, :) = (lactivity - ractivity) ./ (lactivity + ractivity);
    end
    
    weights = 1:num_threshvals;
    numerator = sum(trimmean(TB_LIs, 0.25, 2)' .* weights);
    weighted_li(j,:) = (numerator / sum(weights) * 100);
end
% ft_progress('close');

[~, idx_mx] = max(weighted_li); LI_max = wi(idx_mx,:);

if cfg_main.fplot ==1
    
    figure,plot(weighted_li),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1000   300]);
    %     title([cfg_main.tit, ' - ', tmp.Comment]),
    xlabel('temporal windows (sec)')
    ylabel('LI')
    set(gca,'color','none');
end

end
