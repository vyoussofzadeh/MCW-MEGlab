function [LI,unq_roi_idx] = do_lat_analysis(cfg)

clc, cd(cfg.BS_data_dir), tmp = load(cfg.sinput);

wi = cfg.wi;

idx_L = cfg.lat_index(:,1);
idx_R = cfg.lat_index(:,2);

% thre = 0.5;
LI = []; roi_idx = [];
for j=1:size(wi,1)
    timind1 = nearest(tmp.Time, wi(j,1));
    timind2 = nearest(tmp.Time, wi(j,2));
    [parcelval,rois] = do_sourceparcell_surface(cfg.atlas,mean(tmp.ImageGridAmp(:,timind1:timind2),2));
    %     [roiid, idx, roi_val] = do_barplot_ecp(parcelval,roi_ha2', 0.95, 2);
    [roiid, idx, roi_val] = do_barplot_ecp(parcelval,rois, 0.95, 2);
    %     parcelval(parcelval < thre.*max(parcelval(:))) = 0;
    %     parcelval(parcelval < 0) = nan;
    %     m_left = mean(parcelval(1:2:end)); m_right = mean(parcelval(2:2:end));
    m_left = mean(parcelval(idx_L)); m_right = mean(parcelval(idx_R));
    LI(j) = (m_left - m_right)./ (m_left + m_right);
    roi_idx = [roi_idx, idx];
end

unq_roi_idx = unique(roi_idx);
% % roi_ha2(unq_roi_idx)
% rois(unq_roi_idx)

if cfg.fplot ==1
    figure,plot(LI),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1500   500]);
    title(tmp.Comment), set(gca,'color','none');
end