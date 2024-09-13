function [LI, LI_max, roi_count] = do_lat_analysis_contrast_Counting_parcel(cfg_main)

wi = cfg_main.wi;
atlas = cfg_main.atlas;
idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;
thre = cfg_main.thre;
sinput = cfg_main.sinput;
doavg =  cfg_main.doavg;

% Parcel_based (mean parcels) LI analysis
tmp = load(fullfile(cfg_main.BS_data_dir, sinput));

% removing the negive effects
% tmp.ImageGridAmp(tmp.ImageGridAmp<0) = 0;

%%
% figure, plot(tmp.ImageGridAmp(:,1))

LI = [];
roi_count = [];

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
    [LI_clin, roi_cnt] = do_LI_clincial_parcel(cfg);
    LI(j) = LI_clin;
    roi_count(j) = roi_cnt;
end

if cfg_main.fplot ==1
    figure,plot(LI),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1000   300]);
%     title([cfg_main.tit, ' - ', tmp.Comment]),
    xlabel('temporal windows (sec)')
    ylabel('LI')
    set(gca,'color','none');
end

[~, idx_mx] = max(LI); LI_max = wi(idx_mx,:);

