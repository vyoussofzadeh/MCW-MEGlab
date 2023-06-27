function [mm, wi_sub_max, source_val] = do_corr_analysis(cfg_main)

d_in = cfg_main.d_in;
metric = cfg_main.metric;
var_name = cfg_main.var_name;

%%
metric(isnan(metric(:))) = 0;
% tmp  = S_data_anim_pt_new_tle.d_left.all_metrics; tmp(isnan(tmp(:))) = 0;

var_name_num = [];
for j=1:length(var_name)
    var_name_num{j} = [num2str(j), '-', var_name{j}];
end
var_name_num = strrep(var_name_num, '_', '-');

% size(LI_sub_ltle)
% size(tmp)

% d_in = mean(LI_sub_ltle,2);
% d_in = max(LI_sub_ltle');
size(d_in)

cr_val = [];
for i=1:size(metric,2)
   cr = corrcoef(d_in, metric(:,i));
   cr_val(i) = cr(1,2);
end

figure, bar(cr_val, 0.4)
set(gca,'Xtick', 1:length(cr_val),'XtickLabel',[var_name_num]);
set(gca,'FontSize',8,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   400   1000   500]);
grid

[mx, mxid] = max(cr_val);
figure, plot(d_in,metric(:,mxid),'*'), title(var_name_num(mxid))