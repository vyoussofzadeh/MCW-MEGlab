function [tbl3, idx3, roi_val] = do_barplot_ecp(pow_parcel,rois, thre, plotflag)

L = length(rois);
if plotflag==1
    figure,
    bar(pow_parcel);
    set(gca,'Xtick', 1:L,'XtickLabel',1:L);
    box off
    set(gca,'color','none');
    xlim([0,L+1])
    xlabel('ROI');
    ylabel('Source power');
    set(gcf, 'Position', [1000   100   1200   500]);
    set(gca,'FontSize',10,'XTickLabelRotation',90);
    grid
end

% for i=1:L
%     roiid{i} = [num2str(i),': ',rois{i}];
% end
% disp(roiid')

% [val,idx]  = sort(pow_parcel,'descend');

val = pow_parcel;
val = (val - min(val(:))) ./ (max(val(:)) - min(val(:)));

idx2 = val >= thre.*max(val);
idx3 = idx(idx2);

% roiid(idx(1:nroi))'
% zpow_parcel = zscore(pow_parcel);
% zpow_parcel(idx(1:nroi))'

pow_parcel_thre = zeros(size(pow_parcel));
pow_parcel_thre(idx3) = pow_parcel(idx3);

roi_val = [];
roi_val.parcelval = pow_parcel_thre';
roi_val.rois = rois(idx3)';

tbl1 = table(round(pow_parcel(idx3),2)');
tbl1.Properties.VariableNames = {'t_val'};
tbl2 = cell2table(rois(idx3)');
tbl2.Properties.VariableNames = {'region'};
tbl3 = [tbl1,tbl2];
% disp(tbl3)
