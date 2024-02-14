function [mm, wi_sub_max, source_val] = do_maxLI_sourcemap_diff(cfg_main)

wi = cfg_main.spec.wi;
sid = cfg_main.spec.sid;
data_save_dir = cfg_main.spec.outdir;
s_tag = cfg_main.spec.stag;

%%
mm = []; cfg = []; wi_sub_max = [];

cfg.net_sel_mutiple_label = cfg_main.spec.net_sel_mutiple_label;
cfg.LI_sub = cfg_main.LI_sub;

cfg.wi = wi;
cfg.plotflag = 0;
for i=1:size(cfg.LI_sub,2)
    cfg.subsel = i;
    mm(i,:) = do_plot_sub_LI(cfg);
    [M, midx] = max(mm(i,:));
    wi_sub_max(i,:) = wi(midx,:);
end

%%
s_tag = strrep(s_tag, '_', '-');
% close all

% d_in = mean(mm');
d_in = max(mm') + min(mm');
figure, bar(d_in, 0.4, 'FaceColor',[.5 .5 .5])
L = length(max(mm'));
set(gca,'Xtick', 1:L,'XtickLabel',sid);
set(gca,'FontSize',8,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   400   1000   300]);
set(gca,'color','none');
xlim([0,L+1])
title(s_tag)
grid

for i = 1:numel(d_in)
    textString = sprintf('%.1f', d_in(i));
    text((i), d_in(i), textString, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
end

cfg = []; cfg.outdir = data_save_dir; cfg.filename = ['gmean-', s_tag];  
cfg.type = 'fig'; do_export_fig(cfg)

%%
val = round(mean(wi(:,1),2),2); figure, plot(val, mm'), title(s_tag)
set(gca,'color','none');

cfg = []; cfg.outdir = data_save_dir; cfg.filename = s_tag;  
cfg.type = 'fig'; do_export_fig(cfg)

%%
figure, bar(mean(mm'), 0.4, 'FaceColor',[.5 .5 .5])
L = length(max(mm'));
set(gca,'Xtick', 1:L,'XtickLabel',sid);
set(gca,'FontSize',8,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   400   1000   300]);
set(gca,'color','none');
xlim([0,L+1])
grid

cfg = []; 
cfg.S_data_sel1 = cfg_main.s1;     cfg.S_data_sel2 = cfg_main.s2;
cfg.BS_data_dir = cfg_main.spec.BS_dir;
cfg.wi_sub_max = wi_sub_max; 
cfg.src = cfg_main.spec.src; 
source_val = do_sourcemap_time_sub_optimal_toi_diff(cfg);

cfg = []; cfg.outdir = data_save_dir; cfg.filename = ['smap-', s_tag];  
cfg.type = 'fig'; do_export_fig(cfg)
