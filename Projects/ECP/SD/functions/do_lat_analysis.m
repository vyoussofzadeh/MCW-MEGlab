function [LI,unq_roi_idx, LI_max] = do_lat_analysis(cfg)

wi = cfg.wi;
atlas = cfg.atlas;

if size(cfg.lat_index,2) == 2
    idx_L = cfg.lat_index(:,1);
    idx_R = cfg.lat_index(:,2);
else
    error('check the size of LI indecies')
end
thre = cfg.thre;
sinput = cfg.sinput;

%%
tmp = load(fullfile(cfg.BS_data_dir, sinput));
LI = []; roi_idx = [];
for j=1:size(wi,1)
    
    timind1 = nearest(tmp.Time, wi(j,1)); timind2 = nearest(tmp.Time, wi(j,2));
    [parcelval,rois] = do_sourceparcell_surface(atlas,mean(tmp.ImageGridAmp(:,timind1:timind2),2));
    [~, idx, ~] = do_barplot_ecp(parcelval,rois, thre, 2);
    
    m_left = mean(parcelval(idx_L));
    m_right = mean(parcelval(idx_R));
    
    LI(j) = (m_left - m_right)./ (m_left + m_right);
    roi_idx = [roi_idx, idx];
end
unq_roi_idx = unique(roi_idx);

if cfg.fplot ==1
    figure,plot(LI),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1500   500]);
    title(tmp.Comment), set(gca,'color','none');
end

[~, idx_mx] = max(LI); LI_max = wi(idx_mx,:);
% [~, idx_mn] = min(LI); LI_min = wi(idx_mn,:);