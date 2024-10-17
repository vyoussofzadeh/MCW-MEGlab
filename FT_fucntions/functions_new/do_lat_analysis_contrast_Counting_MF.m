function [LI, LI_max, roi_count] = do_lat_analysis_contrast_Counting_MF(cfg_main)

wi = cfg_main.wi;
atlas = cfg_main.atlas;
idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;
thre = cfg_main.thre;
sinput = cfg_main.sinput;
doavg =  cfg_main.doavg;

%- Parcel_based (mean parcels) LI analysis
tmp = load(fullfile(cfg_main.BS_data_dir, sinput)); 
if  exist('tmp.Value','var'), tmp.ImageGridAmp = tmp.Value; end

% Look up indicies for verticies or (sub)ROIs from HCP atlas
if size(tmp.ImageGridAmp,1) > 360
    sScout = cfg_main.atlas;
    
    % Get left and right subregions from scout data
    LHscout = [];
    for i = 1:length(idx_L)
        LHscout = [LHscout, sScout.Scouts(idx_L(i)).Vertices];
    end
    
    RHscout = [];
    for i = 1:length(idx_R)
        RHscout = [RHscout, sScout.Scouts(idx_R(i)).Vertices];
    end
    idx_LR_updt = [LHscout,RHscout];
else
    idx_LR_updt = [idx_L,idx_R];
end

switch cfg_main.math
    case 'db'
        Fdata = tmp.ImageGridAmp; tidx = tmp.Time < 0; meanBaseline = mean(Fdata(:,tidx),2);
        Fdata = 10 .* log10(abs(bst_bsxfun(@rdivide, Fdata, meanBaseline)));
        tmp.ImageGridAmp = Fdata;
        tmp.ImageGridAmp(tmp.ImageGridAmp < 0) = 0;
    case 'rectif'
        Fdata = tmp.ImageGridAmp;
        Fdata(Fdata < 0) = 0;
        tmp.ImageGridAmp = Fdata;
    case 'rectif_bslnormal'
        Fdata = tmp.ImageGridAmp(idx_LR_updt,:);
        Fdata(Fdata < 0) = 0;
        tidx = tmp.Time < 0; meanBaseline = mean(Fdata(:,tidx),2);
        Fdata = Fdata./meanBaseline;
        tmp.ImageGridAmp(idx_LR_updt,:) = Fdata;
    case 'rectif_zbsl'
        Fdata = tmp.ImageGridAmp(idx_LR_updt,:);
        Fdata(Fdata < 0) = 0;
        tidx = tmp.Time < 0;
        meanBaseline = mean(Fdata(:,tidx),2);
        stdBaseline = std(Fdata(:,tidx)')';
        Fdata = (Fdata - meanBaseline)./stdBaseline;
        tmp.ImageGridAmp(idx_LR_updt,:) = Fdata;
end

%% Global window threshold
mdwin = [];
for j=1:size(wi,1)
    timind1 = nearest(tmp.Time, wi(j,1)); timind2 = nearest(tmp.Time, wi(j,2));
    dwin = tmp.ImageGridAmp(:,timind1:timind2);
    mdwin(j,:) = mean(dwin(idx_LR_updt,:),2);
end
globalwindowthresh = max(mdwin(:));
% figure, plot(wi(:,1),mdwin);

%%
% figure, plot(tmp.ImageGridAmp(:,1))

LI = [];
roi_count_L = []; roi_count_R = [];
roi_count_L_raw = []; roi_count_R_raw = [];

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
    cfg.parcellaion = cfg_main.parcellaion;
    cfg.globalwindowthresh = globalwindowthresh;
    [LI_clin, roi_sum_cnt, roi_cnt] = do_LI_clinical(cfg);
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
    set(gca, 'color', 'none'); % Transparent backgroundend
end

[~, idx_mx] = max(LI); LI_max = wi(idx_mx,:);

roi_count = [];
roi_count.left = roi_count_L;
roi_count.right = roi_count_R;
roi_count.left_raw = roi_count_L_raw;
roi_count.right_raw = roi_count_R_raw;


