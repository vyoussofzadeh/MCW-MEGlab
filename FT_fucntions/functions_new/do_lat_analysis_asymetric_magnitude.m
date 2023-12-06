function [LI, LI_max, pow] = do_lat_analysis_asymetric_magnitude(cfg_main)


wi = cfg_main.wi;
atlas = cfg_main.atlas;
idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;
thre = cfg_main.thre;
sinput = cfg_main.sinput;

%% Parcel_based (mean parcels) LI analysis
tmp = load(fullfile(cfg_main.BS_data_dir, sinput));

%%
% figure, plot(tmp.ImageGridAmp(:,1))

LI = []; 
for j=1:size(wi,1)
    
    timind1 = nearest(tmp.Time, wi(j,1)); timind2 = nearest(tmp.Time, wi(j,2));
    
    cfg = [];
    cfg.thre = thre;
    cfg.atlas = atlas;
    cfg.d_in = mean(tmp.ImageGridAmp(:,timind1:timind2),2);
    cfg.idx_L = idx_L;
    cfg.idx_R = idx_R;
    cfg.Threshtype = cfg_main.Threshtype;
    cfg.thre = thre;
    [LI_clin,pow] = do_LI_magnitude(cfg);
    LI(j) = LI_clin;
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















% wi = cfg_main.wi;
% atlas = cfg_main.atlas;
% idx_L = cfg_main.index_L;
% idx_R = cfg_main.index_R;
% 
% thre = cfg_main.thre;
% sinput = cfg_main.sinput;
% 
% tmp = load(fullfile(cfg_main.BS_data_dir, sinput));
% LI = zeros(size(wi,1), 1);
% pow_left = zeros(size(wi,1), 1);
% pow_right = zeros(size(wi,1), 1);
% 
% for j = 1:size(wi,1)
%     timind1 = nearest(tmp.Time, wi(j,1));
%     timind2 = nearest(tmp.Time, wi(j,2));
% 
%     cfg = struct('thre', thre, 'atlas', atlas, 'd_in', mean(tmp.ImageGridAmp(:,timind1:timind2),2));
%     [parcelval, ~] = do_sourceparcell_surface(cfg);
% 
%     cfg = struct('thre', thre, 'do_atan', 1, 'do_normal', 0, 'do_plot', 0, 'idx_L', idx_L, 'idx_R', idx_R);
%     parcelthre = do_apply_thre(cfg, parcelval);
% 
%     m_left = nanmean(parcelthre(idx_L));
%     m_right = nanmean(parcelthre(idx_R));
% 
%     LI(j) = 100 * (m_left - m_right) / (m_left + m_right);
%     pow_left(j) = m_left;
%     pow_right(j) = m_right;
% end
% 
% if cfg_main.fplot == 1
%     figure, plot(LI);
%     val = round(mean(wi(:,1), 2), 2);
%     set(gca, 'Xtick', 1:2:length(wi), 'XtickLabel', val(1:2:end));
%     set(gca, 'FontSize', 8, 'XTickLabelRotation', 90);
%     set(gcf, 'Position', [1000, 400, 1000, 300]);
%     title([cfg_main.tit, ' - ', tmp.Comment]);
%     xlabel('temporal windows (sec)');
%     ylabel('LI');
%     set(gca, 'color', 'none');
% end
% 
% [~, idx_mx] = max(LI);
% LI_max = wi(idx_mx, :);
% 
% pow = struct('left', pow_left, 'right', pow_right);
% 
% end
