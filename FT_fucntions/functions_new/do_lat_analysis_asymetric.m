function [LI, LI_max] = do_lat_analysis_asymetric(cfg_main)

wi = cfg_main.wi;
atlas = cfg_main.atlas;
idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;

thre = cfg_main.thre;
sinput = cfg_main.sinput;

%% Parcel_based (mean parcels) LI analysis
tmp = load(fullfile(cfg_main.BS_data_dir, sinput));
LI = []; roi_idx = [];
for j=1:size(wi,1)

    timind1 = nearest(tmp.Time, wi(j,1)); timind2 = nearest(tmp.Time, wi(j,2));
    [parcelval, ~] = do_sourceparcell_surface(atlas,mean(tmp.ImageGridAmp(:,timind1:timind2),2));
    
    cfg = []; 
    cfg.thre = thre; 
    cfg.do_atan = 1; 
    cfg.do_normal = 0; 
    cfg.do_plot = 0;
    cfg.idx_L = idx_L; % region-based threshold
    cfg.idx_R = idx_R; % region-based threshold
    parcelthre = do_apply_thre(cfg, parcelval);
    
%     [~, idx, parcelthre] = do_barplot_ecp(parcelval,rois, thre, 2);
%     parcelval = parcelthre.parcelval;

    m_left = nanmean(parcelthre(idx_L));
    m_right = nanmean(parcelthre(idx_R));

    LI(j) = 100*(m_left - m_right)./ (m_left + m_right);
%     roi_idx = [roi_idx, idx];
end
% unq_roi_idx = unique(roi_idx);

if cfg_main.fplot ==1
    figure,plot(LI),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1000   300]);
    title([cfg_main.tit, ' - ', tmp.Comment]),
    xlabel('temporal windows (sec)')
    ylabel('LI')
    set(gca,'color','none');
end

[~, idx_mx] = max(LI); LI_max = wi(idx_mx,:);
% [~, idx_mn] = min(LI); LI_min = wi(idx_mn,:);

%% vertex LI analysis
% tmp = load(fullfile(cfg.BS_data_dir, sinput));
% LI = [];
% % roi_idx = [];
% for j=1:size(wi,1)
%     
%     timind1 = nearest(tmp.Time, wi(j,1)); timind2 = nearest(tmp.Time, wi(j,2));
%     tmp2 = mean(tmp.ImageGridAmp(:,timind1:timind2),2);
%     %     tmp2(tmp2 < thre.*max(tmp2(:))) = 0;
%     [parcelval, pow_parcel_max, rois] = do_source_surface(atlas,tmp2);
% 
%     parcelvalmat_L = cat(1,parcelval{idx_L});
%     parcelvalmat_R = cat(1,parcelval{idx_R});
%     
%     mx = max([max(parcelvalmat_L), max(parcelvalmat_R)]);
%     
% %     figure, plot(parcelvalmat_L), hold on, yline(thre.*mx)
% %     figure, plot(parcelvalmat_R), hold on, yline(thre.*mx)
%     
%     parcelvalmat_L = parcelvalmat_L(parcelvalmat_L > thre.*mx);
%     parcelvalmat_R = parcelvalmat_R(parcelvalmat_R > thre.*mx);
%     
%     m_left = nanmean(parcelvalmat_L(idx_L));
%     m_right = nanmean(parcelvalmat_R(idx_R));
%     
%     LI(j) = 100*(m_left - m_right)./ (m_left + m_right);
% 
%     %%
% %     LI(j) = 100*(length(parcelvalmat_L) - length(parcelvalmat_R))/ (length(parcelvalmat_L) + length(parcelvalmat_R));
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     %%
%     
% %     threshold_global = thre.*max(cell2mat(pow_parcel_max));
% %     threshold_global = 0;
% %     %     parcelval{1}
% %     %     pause
% %     
% %     m_left = 0; m_right = 0;
% %     for i=1:length(idx_L)
% %         LHvals  = parcelval{idx_L(i)};
% % %         figure, plot(LHvals), 
% %         threshold  =  thre.*max(LHvals); %pow_parcel_max{idx_L(i)};
% % %         hold on, yline(threshold_global)
% %         ind_L = find(LHvals > threshold_global);
% %         L_ROIcount = length(ind_L);
% %         m_left = m_left + L_ROIcount;
% %         
% %         RHvals  = parcelval{idx_R(i)};
% % %         threshold  =  thre.*max(RHvals); %pow_parcel_max{idx_L(i)};
% %         ind_R = find(RHvals > threshold_global);
% %         R_ROIcount = length(ind_R);
% %         m_right = m_right + R_ROIcount;
% %     end
% %     
% % %     for i=1:length(idx_R)
% % %         parcelval_sel  = parcelval{idx_R(i)};
% % %         parcelval_sel_max  =  max(parcelval_sel); %pow_parcel_max{idx_L(i)};
% % %         m_right = m_right + length(find(parcelval_sel > thre.*parcelval_sel_max));
% % %     end
% %     
% %     %     [~, idx, parcelthre] = do_barplot_ecp(parcelval,rois, thre, 2);
% %     %     parcelval = parcelthre.parcelval;
% %     %
% %     %     m_left = nanmean(parcelval(idx_L));
% %     %     m_right = nanmean(parcelval(idx_R));
% %     
% %     %     m_left = length(parcelval(idx_L) > 0);
% %     %     m_right = length(parcelval(idx_R) > 0);
% %     
% %     LI(j) = 100*(m_left - m_right)./ (m_left + m_right);
%     %     roi_idx = [roi_idx, idx];
% end
% % unq_roi_idx = unique(roi_idx);
% 
% if cfg.fplot ==1
%     figure,plot(LI),
%     val = round(mean(wi(:,1),2),2);
%     set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
%     set(gca,'FontSize',8,'XTickLabelRotation',90);
%     set(gcf, 'Position', [1000   400   1000   300]);
%     title([cfg.tit, ' - ', tmp.Comment]),
%     xlabel('temporal windows (sec)')
%     ylabel('LI')
%     set(gca,'color','none');
% end
% 
% [~, idx_mx] = max(LI); LI_max = wi(idx_mx,:);
% [~, idx_mn] = min(LI); LI_min = wi(idx_mn,:);