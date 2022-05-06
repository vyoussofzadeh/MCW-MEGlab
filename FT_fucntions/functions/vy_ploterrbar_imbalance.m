function vy_ploterrbar_imbalance(cfg, datain1, datain2)

bar([mean(datain1); mean(datain2)]', 0.4), title(cfg.titletag)
hold on
errorbar(1:2,[mean(datain1); mean(datain2)]',[std(datain1); std(datain2)]')
set(gca,'Xtick', 1:2,'XtickLabel',cfg.XtickLabeltag);
xlabel(cfg.xlabeltag), ylabel(cfg.ylabeltag),
% set(gca,'FontSize',10,'XTickLabelRotation',90);
grid on
set(gca,'color','none');
% legend({'Str';'Math'})

disp('==============')
disp([num2str(mean(datain1)),' +- ', num2str(std(datain1))])
disp([num2str(mean(datain2)),' +- ', num2str(std(datain2))])
disp([num2str(mean([datain1, datain2])),' +- ', num2str(std([datain1, datain2]))])
disp('==============')
