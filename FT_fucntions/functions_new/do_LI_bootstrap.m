function [num_threshvals, weighted_li] = do_LI_bootstrap(cfg_main)

wi = cfg_main.wi;
atlas = cfg_main.atlas;
idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;

% thre = cfg_main.thre;
sinput = cfg_main.sinput;

RESAMPLE_RATIO = 0.75;
MIN_NUM_THRESH_VOXELS = round(5 / RESAMPLE_RATIO);

tmp = load(fullfile(cfg_main.BS_data_dir, sinput));
sScout = atlas;

% weighted_li = zeros(size(wi, 1), 1);
weighted_li = [];
% weighted_li = zeros(size(wi, 1), size(threshvals, 2));

for j = 1:size(wi, 1)
    
    timind1 = nearest(tmp.Time, wi(j,1));
    timind2 = nearest(tmp.Time, wi(j,2));
    
    d_in = mean(tmp.ImageGridAmp(:, timind1:timind2), 2);
    ImageGridAmp = abs(d_in);
    
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
    LHvals = ImageGridAmp(LHscout);
    RHvals = ImageGridAmp(RHscout);
    
    % Calculate maximum values for left and right subregions
    LH_max = max(LHvals(:));
    RH_max = max(RHvals(:));
    ROIMax = max(LH_max, RH_max);
    
    divs = 2;
    lvals_nonnegative = LHvals(LHvals >= 0);
    rvals_nonnegative = RHvals(RHvals >= 0);
    
    threshvals = (0:(divs-1)) * (ROIMax / (divs - 1));
    
    lr_threshvals = cell(1, numel(threshvals));
    num_threshvals = 0;
    
    for i = 1:numel(threshvals)
        threshval = threshvals(i);
        l_threshvals = lvals_nonnegative(lvals_nonnegative >= threshval);
        r_threshvals = rvals_nonnegative(rvals_nonnegative >= threshval);
        
        if (numel(l_threshvals) >= MIN_NUM_THRESH_VOXELS && numel(r_threshvals) >= MIN_NUM_THRESH_VOXELS)
            num_threshvals = num_threshvals + 1;
            lr_threshvals{num_threshvals} = {l_threshvals, r_threshvals};
        else
            break;
        end
    end
    
    reps_per_li = 5 * num_threshvals;
    TB_LIs = zeros(num_threshvals, reps_per_li);
    
    for i = 1:num_threshvals
        l_set = lr_threshvals{i}{1};
        r_set = lr_threshvals{i}{2};
        
        l_n = max(round(numel(l_set) * RESAMPLE_RATIO), MIN_NUM_THRESH_VOXELS);
        r_n = max(round(numel(r_set) * RESAMPLE_RATIO), MIN_NUM_THRESH_VOXELS);
        
        l_indices = randi(numel(l_set), reps_per_li, l_n);
        lactivity = sum(l_set(l_indices), 2);
        
        r_indices = randi(numel(r_set), reps_per_li, r_n);
        ractivity = sum(r_set(r_indices), 2);
        
        TB_LIs(i, :) = (lactivity - ractivity) ./ (lactivity + ractivity);
    end
    
    weights = 1:num_threshvals;
    numerator = sum(trimmean(TB_LIs, 0.25, 2) .* weights);
    weighted_li(j,:) = (numerator / sum(weights) * 100);
    %     disp(j)
end

if cfg_main.fplot ==1
    
    figure,plot(weighted_li),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1000   300]);
    title([cfg_main.tit, ' - ', tmp.Comment]),
    xlabel('temporal windows (sec)')
    ylabel('LI')
    set(gca,'color','none');
end

end
