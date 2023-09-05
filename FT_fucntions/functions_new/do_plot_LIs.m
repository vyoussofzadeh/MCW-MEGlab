function [wi]  = do_plot_LIs(cfg_main)

DataArray = cfg_main.DataArray;
tit = cfg_main.title;
outdir = cfg_main.outdir;

% DataArray = [megLI_sub_pt, fmri_LIs_val];
Colors = [0.4922    0.0039    0.9063;0.9922    0.8672    0.0039];
figure
% subplot 121
[xPositions, yPositions, ~, ~] = UnivarScatter(DataArray,'Label',{'meg','fMRI'},'MarkerFaceColor',Colors);
ylabel('Laterality','FontSize', 16);
xlabel('Modality','FontSize', 16);
set(gca,'color','none');
title(tit,'fontsize',16)
disp(['fmri: ',num2str(mean(DataArray(:,1))), '+-', num2str(std(DataArray(:,1)))]);
disp(['meg: ',num2str(mean(DataArray(:,2))), '+-', num2str(std(DataArray(:,2)))]);
set(gca,'FontName','HelveticaNeueLT Std Lt');
hold on
f = [xPositions, yPositions];
for j=1:length(f)
    line([f(j,1),f(j,2)],[f(j,3),f(j,4)]);
end

if cfg_main.savefig == 1
    % - export figs
    cfg = []; cfg.outdir = outdir; cfg.filename = tit;
    cfg.type = 'fig'; do_export_fig(cfg)
end

end