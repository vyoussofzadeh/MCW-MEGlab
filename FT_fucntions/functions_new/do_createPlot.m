% New Helper function for plotting
function do_createPlot(yData, val, color, net_sel_mutiple_label, tag, yLabel)

plot(yData,'LineWidth',3, 'color', color),
set(gca,'Xtick', 1:2:length(val),'XtickLabel',val(1:2:end));
set(gca,'FontSize',8,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   400   1100   500]);
lgnd = legend([net_sel_mutiple_label; 'mean']);
title(tag)
ylabel(yLabel)
xlabel('time')
set(gca,'color','none');
set(lgnd,'color','none');

end

