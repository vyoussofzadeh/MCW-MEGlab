function [LI, LI_max, roi_count] = do_lat_analysis_baseline_Counting_MF(cfg_main)

wi = cfg_main.wi;
atlas = cfg_main.atlas;
idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;
thre = cfg_main.thre;
sinput = cfg_main.sinput;
% Threshtype = cfg_main.Threshtype;
doavg = cfg_main.doavg;

%% Parcel_based (mean parcels) LI analysis
tmp_1 = load(fullfile(cfg_main.BS_data_dir, sinput{1})); tmp_1.ImageGridAmp = tmp_1.Value;
tmp_2 = load(fullfile(cfg_main.BS_data_dir, sinput{2})); tmp_2.ImageGridAmp = tmp_2.Value;

idx_LR = [idx_L,idx_R];


% disp({tmp_1.Comment; tmp_2.Comment})

tmp = tmp_1;

switch cfg_main.math
    case 'db'
        tmp.ImageGridAmp = 10 * log10(tmp_1.ImageGridAmp) - 10 * log10(tmp_2.ImageGridAmp); tmp.ImageGridAmp(tmp.ImageGridAmp < 0) = 0;
        
    case 'db_baseline'
        Fdata1 = tmp_1.ImageGridAmp; tidx = tmp_1.Time < 0; meanBaseline = mean(Fdata1(:,tidx),2);
        Fdata1 = 10 .* log10(abs(bst_bsxfun(@rdivide, Fdata1, meanBaseline)));
        
        Fdata2 = tmp_2.ImageGridAmp; tidx = tmp_2.Time < 0; meanBaseline = mean(Fdata2(:,tidx),2);
        Fdata2 = 10 .* log10(abs(bst_bsxfun(@rdivide, Fdata2, meanBaseline)));
        
        tmp.ImageGridAmp = Fdata1 - Fdata2;
        tmp.ImageGridAmp(tmp.ImageGridAmp < 0) = 0;
        
    case 'diff'
        tmp.ImageGridAmp = tmp_1.ImageGridAmp - tmp_2.ImageGridAmp;
        
    case 'rec_diff'
        tmp.ImageGridAmp = tmp_1.ImageGridAmp - tmp_2.ImageGridAmp;
        tmp.ImageGridAmp(tmp.ImageGridAmp < 0) = 0;
        %         tmp.ImageGridAmp = tmp.ImageGridAmp + 0.1;
        
    case 'rec_diff_mbsl'
        tmp.ImageGridAmp(idx_LR,:) = tmp_1.ImageGridAmp(idx_LR,:) - tmp_2.ImageGridAmp(idx_LR,:);
        Fdata = tmp.ImageGridAmp(idx_LR,:); Fdata(Fdata < 0) = 0;
        tidx = tmp.Time < 0; meanBaseline = mean(Fdata(:,tidx),2);
        Fdata = Fdata./meanBaseline; tmp.ImageGridAmp(idx_LR,:) = Fdata;
    case 'bsl_diff_rec'
        Fdata = tmp_1.ImageGridAmp(idx_LR,:);
        tidx = tmp.Time < 0; meanBaseline = mean(Fdata(:,tidx),2); Fdata = Fdata./meanBaseline; % bsl norm
        tmp_1.ImageGridAmp(idx_LR,:) = Fdata;
        
        Fdata = tmp_2.ImageGridAmp(idx_LR,:);
        tidx = tmp.Time < 0; meanBaseline = mean(Fdata(:,tidx),2); Fdata = Fdata./meanBaseline; % bsl norm
        tmp_2.ImageGridAmp(idx_LR,:) = Fdata;
        
        tmp.ImageGridAmp(idx_LR,:) = tmp_1.ImageGridAmp(idx_LR,:) - tmp_2.ImageGridAmp(idx_LR,:);
        Fdata = tmp.ImageGridAmp(idx_LR,:);
        Fdata(Fdata < 0) = 0;
        tmp.ImageGridAmp(idx_LR,:) = Fdata;
end

LI = [];
roi_count_L = []; roi_count_R = [];


for j=1:size(wi,1)
    
    timind1 = nearest(tmp.Time, wi(j,1)); timind2 = nearest(tmp.Time, wi(j,2));
    
    cfg = [];
    cfg.thre = thre;
    cfg.atlas = atlas;
    if doavg == 1
        cfg.d_in = mean(tmp.ImageGridAmp(:,timind1:timind2),2);
    else
        cfg.d_in = tmp.ImageGridAmp(:,timind1:timind2);
    end
    cfg.idx_L = idx_L;
    cfg.idx_R = idx_R;
    cfg.Threshtype = cfg_main.Threshtype;
    cfg.thre = thre;
    cfg.parcellaion = cfg_main.parcellaion;
    [LI_clin, roi_sum_cnt, roi_cnt] = do_LI_clincial(cfg);
    LI(j) = LI_clin;
    roi_count_sum(j) = roi_sum_cnt;
    roi_count_L(j) = roi_cnt.L;
    roi_count_R(j) = roi_cnt.R;
end

if cfg_main.fplot ==1
    figure;
    yyaxis left; % Left y-axis for LI
    plot(LI,'LineWidth',1.5);
    ylabel('Lateralization Index (LI)');
    
    hold on;
    
    yyaxis right; % Right y-axis for Power
    plot(roi_count_L,'LineWidth',1.5); % Mean power for left activities
    plot(roi_count_R,'LineWidth',1.5); % Mean power for right activities
    ylabel('Count');
    
    legend({'LI', 'Left Count', 'Right Count'}, 'Location', 'best');
    
    % Set x-axis ticks and labels
    val = round(mean(wi(:,1),2),2);
    set(gca, 'Xtick', 1:2:length(wi), 'XtickLabel', val(1:2:end));
    set(gca, 'FontSize', 8, 'XTickLabelRotation', 90);
    set(gcf, 'Position', [1000, 400, 1000, 300]);
    
    xlabel('Mean Temporal Windows (sec)');
    title(['Counting']);
    set(gca, 'color', 'none'); % Transparent background
end

[~, idx_mx] = max(LI); LI_max = wi(idx_mx,:);

roi_count = [];
roi_count.left = roi_count_L;
roi_count.right = roi_count_R;

