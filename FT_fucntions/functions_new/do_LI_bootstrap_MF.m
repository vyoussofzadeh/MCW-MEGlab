function [weighted_li, LI_max, num_valid_thresh, roi_count] = do_LI_bootstrap_MF(cfg_main)

% MIN_NUM_THRESH_VOXELS is a predefined minimum number ...
% of voxels that are required for the computation to be considered valid.

wi = cfg_main.wi;
atlas = cfg_main.atlas;
idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;
divs = cfg_main.divs; % divs: number of divisions between 0 and max value in img
n_resampling = cfg_main.n_resampling; % 200
doavg = cfg_main.doavg;
downsamplerate = cfg_main.downsamplerate;
atlas = cfg_main.atlas;


sinput = cfg_main.sinput;

RESAMPLE_RATIO = cfg_main.RESAMPLE_RATIO; % 0.75
MIN_NUM_THRESH_VOXELS = round(5 / RESAMPLE_RATIO);

%%
tmp = load(fullfile(cfg_main.BS_data_dir, sinput)); tmp.ImageGridAmp = tmp.Value;

idx_LR = [idx_L,idx_R];

switch cfg_main.math
    case 'db'
        Fdata = tmp.ImageGridAmp; tidx = tmp.Time < 0; meanBaseline = mean(Fdata(:,tidx),2);
        Fdata = 10 .* log10(abs(bst_bsxfun(@rdivide, Fdata, meanBaseline)));
        tmp.ImageGridAmp = Fdata;
        tmp.ImageGridAmp(tmp.ImageGridAmp < 0) = 0;
    case 'raw'
        %         tmp.ImageGridAmp(tmp.ImageGridAmp < 0) = 0;
    case 'rectif'
        Fdata = tmp.ImageGridAmp;
        Fdata(Fdata < 0) = 0;
        %         Fdata = Fdata + 0.1;
        tmp.ImageGridAmp = Fdata;
    case 'rectif_bslnormal'
        Fdata = tmp.ImageGridAmp(idx_LR,:);
        Fdata(Fdata < 0) = 0;
        tidx = tmp.Time < 0; meanBaseline = mean(Fdata(:,tidx),2);
        Fdata = Fdata./meanBaseline;
        tmp.ImageGridAmp(idx_LR,:) = Fdata;
    case 'rectif_bslnormal_mean'
        
        Fdata_left  = mean(tmp.ImageGridAmp(idx_L,:),1);
        Fdata_right = mean(tmp.ImageGridAmp(idx_R,:),1);
        
        tidx = tmp.Time < 0;
        
        mbsl_L = mean(Fdata_left(:,tidx),2); Fdata_left = Fdata_left./mbsl_L;
        mbsl_R = mean(Fdata_right(:,tidx),2); Fdata_right = Fdata_right./mbsl_R;
        
        for j = 1:length(idx_L), tmp.ImageGridAmp(idx_L(j),:) = Fdata_left; end
        for j = 1:length(idx_R), tmp.ImageGridAmp(idx_R(j),:) = Fdata_right; end
end

%%
sScout = atlas;

weighted_li = [];
% ft_progress('init', 'text',     'please wait ...');

% Predefine the arrays with max size, will truncate later
left_activities = zeros(size(wi, 1), divs);
right_activities = zeros(size(wi, 1), divs);

for j = 1:size(wi, 1)
    
    %     ft_progress(j/size(wi, 1), 'Processing interval %d from %d', j, size(wi, 1));
    timind1 = nearest(tmp.Time, wi(j,1));
    timind2 = nearest(tmp.Time, wi(j,2));
    
    if doavg ==  1
        d_in = mean(tmp.ImageGridAmp(:, timind1:timind2), 2);
    else
        d_in = tmp.ImageGridAmp(:, timind1:timind2);
    end
    
    ImageGridAmp = abs(d_in);
    
    % downsample data
    ImageGridAmp = downsample(ImageGridAmp',downsamplerate);
    ImageGridAmp = ImageGridAmp';
    
    if size(ImageGridAmp,1) > 360
        cfg_main.parcellaion = 0;
    end
    
    if cfg_main.parcellaion == 1
        % Extract amplitude values for left and right subregions
        LHvals = ImageGridAmp(idx_L,:);
        RHvals = ImageGridAmp(idx_R,:);
        
    else
        % Get left and right subregions from scout data
        LHscout = [];
        for i = 1:length(idx_L)
            LHscout = [LHscout, sScout.Scouts(idx_L(i)).Vertices];
        end
        
        RHscout = [];
        for i = 1:length(idx_R)
            RHscout = [RHscout, sScout.Scouts(idx_R(i)).Vertices];
        end
        
        % Extract amplitude values for left and right subregions
        LHvals = ImageGridAmp(LHscout,:);
        RHvals = ImageGridAmp(RHscout,:);
    end
    
    % Calculate maximum values for left and right subregions
    LH_max = max(LHvals(:));
    RH_max = max(RHvals(:));
    ROIMax = max(LH_max, RH_max);
    
    threshvals = (0:(divs-1)) * (ROIMax / (divs - 1));
    l_threshvals_all = LHvals(LHvals(:) >= 0);
    r_threshvals_all = RHvals(RHvals(:) >= 0);
    
    num_l_valid_thresh = sum(arrayfun(@(x) sum(l_threshvals_all >= x), threshvals) >= MIN_NUM_THRESH_VOXELS);
    num_r_valid_thresh = sum(arrayfun(@(x) sum(r_threshvals_all >= x), threshvals) >= MIN_NUM_THRESH_VOXELS);
    
    % Get the minimum count of valid thresholds for both hemispheres for the current window
    num_valid_thresh = min(num_l_valid_thresh, num_r_valid_thresh);
    valid_threshvals = threshvals(1:num_valid_thresh);
    
    % Initialize local arrays to store activities for current window
    left_act = zeros(1, num_valid_thresh);
    right_act = zeros(1, num_valid_thresh);
    
    lr_threshvals = cell(1, num_valid_thresh);
    
    for i = 1:num_valid_thresh
        threshval = valid_threshvals(i);
        l_threshvals = l_threshvals_all(l_threshvals_all >= threshval);
        r_threshvals = r_threshvals_all(r_threshvals_all >= threshval);
        lr_threshvals{i} = {l_threshvals, r_threshvals};
    end
    
    TB_LIs = zeros(num_valid_thresh, n_resampling);
    
    for i = 1:num_valid_thresh
        l_set = lr_threshvals{i}{1};
        r_set = lr_threshvals{i}{2};
        
        l_n = max(round(numel(l_set) * RESAMPLE_RATIO), MIN_NUM_THRESH_VOXELS);
        r_n = max(round(numel(r_set) * RESAMPLE_RATIO), MIN_NUM_THRESH_VOXELS);
        
        l_indices = randi(numel(l_set), n_resampling, l_n);
        lactivity = sum(l_set(l_indices), 2);
        
        r_indices = randi(numel(r_set), n_resampling, r_n);
        ractivity = sum(r_set(r_indices), 2);
        
        % Store summed activities
        left_act(i) = mean(lactivity);
        right_act(i) = mean(ractivity);
        
        TB_LIs(i, :) = (lactivity - ractivity) ./ (lactivity + ractivity);
    end
    
    % Store the activities for each window
    left_activities(j, 1:num_valid_thresh) = left_act;
    right_activities(j, 1:num_valid_thresh) = right_act;
    
    weights = 1:num_valid_thresh;
    numerator = sum(trimmean(TB_LIs, 0.25, 2)' .* weights);
    weighted_li(j,:) = (numerator / sum(weights) * 100);
end
% ft_progress('close');

roi_count = [];
roi_count.left = left_activities;
roi_count.right = right_activities;

[~, idx_mx] = max(weighted_li); LI_max = wi(idx_mx,:);

if cfg_main.fplot ==1
    
    figure;
    yyaxis left; % Left y-axis for LI
    plot(weighted_li,'LineWidth',1.5);
    ylabel('Lateralization Index (LI)');
    
    hold on;
    
    yyaxis right; % Right y-axis for Power
    plot(mean(left_activities,2),'LineWidth',1.5); % Mean power for left activities
    plot(mean(right_activities,2),'LineWidth',1.5); % Mean power for right activities
    ylabel('Mean Count');
    
    legend({'LI', 'Left Count', 'Right Count'}, 'Location', 'best');
    
    % Set x-axis ticks and labels
    val = round(mean(wi(:,1),2),2);
    set(gca, 'Xtick', 1:2:length(wi), 'XtickLabel', val(1:2:end));
    set(gca, 'FontSize', 8, 'XTickLabelRotation', 90);
    set(gcf, 'Position', [1000, 400, 1000, 300]);
    
    xlabel('Mean Temporal Windows (sec)');
    title(['Bootstrap']);
    set(gca, 'color', 'none'); % Transparent background
    
end

end

