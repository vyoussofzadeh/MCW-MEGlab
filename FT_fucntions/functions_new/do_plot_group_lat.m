
% Additional functions for plotLINetworks, plotMeanLIWithStdDev, plotGeneralMeanLI, and saveFigure
% can be defined similarly with the relevant code blocks from the original script.


function do_plot_group_lat(cfg_main)


LI_sub = cfg_main.LI_sub;
wi = cfg_main.wi;
net_sel_mutiple_label = cfg_main.net_sel_mutiple_label;
S_data_sel = cfg_main.S_data_sel;
outdir = cfg_main.outdir;
network_sel = cfg_main.network_sel;

%%
colr = distinguishable_colors(length(network_sel));

mLI_sub = squeeze(nanmean(LI_sub,2));
figure,
clear LI_val
for j=1:length(network_sel)
    LI_val(j,:) =  mLI_sub(network_sel(j),:);
    plot(LI_val(j,:), 'Color', colr(j,:)),
    val = round(nanmean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   700   800   300]);
    hold on

%     [peak_y, peak_x] = max(LI_val(j,:));
%     xline(peak_x, 'color', colr(j,:),'LineWidth',1);
%     val = round(mean(wi(:,1),2),2);
%     x = val(peak_x);
%     txt = num2str(x);
%     text(peak_x, peak_y, txt, 'Rotation', 90)
end
plot(nanmean(LI_val(:,:)),'LineWidth',3);
lgnd = legend([net_sel_mutiple_label(network_sel); 'mean']);
legend('AutoUpdate', 'off')
set(lgnd,'color','none');
% title(['mean LIs: ', S_data_sel.s_tag])
set(lgnd,'color','none');
set(gca,'color','none');
xlabel('time')
ylabel('LI')

for j=1:length(network_sel)
    LI_val(j,:) =  mLI_sub(network_sel(j),:);
    [peak_y, peak_x] = max(LI_val(j,:));
    xline(peak_x, 'color', colr(j,:),'LineWidth',1);
    val = round(mean(wi(:,1),2),2);
    x = val(peak_x);
    txt = num2str(x);
    text(peak_x, peak_y, txt, 'Rotation', 90)
end

if cfg_main.savefig == 1
    cfg = [];
    cfg.outdir = outdir;
    filename = ['meanLIs_', S_data_sel.s_tag];
    cfg.filename = filename;
    cfg.type = 'fig';
    do_export_fig(cfg)
end

%%
nn = round(length(network_sel)/3);

clear std_dev
figure,
for j=1:length(network_sel)
    tmp = squeeze(LI_sub(network_sel(j),:,:));
%     subplot(nn,3,j)
    plot(tmp'),
    hold on
    mm(j,:) = nanmean(tmp);
    std_dev(j,:) = std(tmp);
    plot(mm(j,:),'LineWidth',3, 'color', colr(j,:)),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:5:length(wi),'XtickLabel',val(1:5:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gca,'color','none');
%     title([S_data_sel.s_tag, '-', net_sel_mutiple_label{network_sel(j)}])
    ylim([-110, 110])
    [peak_y, peak_x] = max(mm(j,:));
    xline(peak_x, 'color', colr(j,:),'LineWidth',1);
    val = round(mean(wi(:,1),2),2);
    x = val(peak_x);
    txt = num2str(x);
    text(peak_x, peak_y, txt, 'Rotation', 90)
end

set(gcf, 'Position', [1000   700   800   300]);
if cfg_main.savefig == 1
    cfg = [];
    cfg.outdir = outdir;
    filename = ['meanLIs_ROIs_', S_data_sel.s_tag];
    cfg.filename = filename;
    cfg.type = 'fig';
    do_export_fig(cfg)
end

%% plot mean LI
% figure,
% for j=1:size(mm,1)
%     hold on
%     plot(mm(j,:),'LineWidth',1, 'color', colr(j,:)),
%     val = round(mean(wi(:,1),2),2);
%     set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
%     set(gca,'FontSize',8,'XTickLabelRotation',90);
%     set(gcf, 'Position', [1000   400   1100   500]);
% end
% lgnd = legend(net_sel_mutiple_label(network_sel));
% title(['subject LIs: ', S_data_sel.s_tag,])
% set(gca,'color','none');
% set(lgnd,'color','none');
% legend('AutoUpdate', 'off')
% for j=1:size(mm,1)
%     hold on
%     [peak_y, peak_x] = max(mm(j,:));
%     xline(peak_x, 'color', colr(j,:),'LineWidth',1);
%     val = round(mean(wi(:,1),2),2);
%     x = val(peak_x);
%     txt = num2str(x);
% %     text(peak_x, peak_y, txt, 'Rotation', 90)
% end
%
% if cfg_main.savefig == 1
%     cfg = [];
%     cfg.outdir = outdir;
%     filename = ['meanLIs_meanROIs_', S_data_sel.s_tag];
%     cfg.filename = filename;
%     cfg.type = 'fig';
%     do_export_fig(cfg)
% end

%%
% close all
figure,
for j=1:size(mm,1)
    a = mm(j,:);
    b = 1:numel(a);
    curve1 = a + std_dev(j,:);
    curve2 = a - std_dev(j,:);

    x2 = [b, fliplr(b)];
    inBetween = [curve1, fliplr(curve2)];
    hold on;
    h = plot(b, a, 'color', colr(j,:), 'LineWidth', 2);
    h2 = fill(x2, inBetween, colr(j,:), 'EdgeColor', 'none', 'facealpha', 0.1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off'; % make the legend for step plot off

end
lgnd = legend(net_sel_mutiple_label(network_sel));
set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
set(gca,'FontSize',8,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   700   800   300]);
% title(['subject LIs: ', S_data_sel.s_tag,])
set(gca,'color','none');
set(lgnd,'color','none');
ylabel('LI')
xlabel('Time')

if cfg_main.savefig == 1
    cfg = [];
    cfg.outdir = outdir;
    filename = ['meanLIs_meanROIs_withvar_', S_data_sel.s_tag];
    cfg.filename = filename;
    cfg.type = 'fig';
    do_export_fig(cfg)
end

%%
d_in = mean(mean(LI_sub,1),3); L = length(d_in);
figure, bar(d_in,0.4)
set(gca,'Xtick', 1:length(d_in),'XtickLabel',1:length(d_in));
set(gca,'FontSize',8,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   700   800   300]);
set(gca,'color','none');
title('mean LI')
xlabel('Subj')
set(lgnd,'color','none');
ylabel('Laterality')
grid

if cfg_main.savefig == 1
    cfg = [];
    cfg.outdir = outdir;
    filename = ['gmeanLIs_', S_data_sel.s_tag];
    cfg.filename = filename;
    cfg.type = 'fig';
    do_export_fig(cfg)
end

end


