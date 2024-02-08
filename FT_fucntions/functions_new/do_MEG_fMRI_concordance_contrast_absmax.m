function [megLI_sub_pt, fmri_LIs_val] = do_MEG_fMRI_concordance_contrast_absmax(cfg_main)


wi = cfg_main.wi;
LI_pt_val_new = cfg_main.LI_val;
net_sel = cfg_main.net_sel;
fmri_LIs_val = cfg_main.fmri_LIs_val;
% net_sel_mutiple_label = cfg_main.net_sel_mutiple_label;
ID = cfg_main.ID;
thre = cfg_main.thre;
savefig = cfg_main.savefig;
outdir = cfg_main.outdir;
absmax_idx = cfg_main.absmax;
bf = cfg_main.buffervalue;

% max_time_idx = absmax_idx.max_time_idx;
max_time_idx = absmax_idx.max_time_pts_sel;


mLI_sub1 = [];
for j=1:length(max_time_idx)
    mLI_sub1(j,:) = LI_pt_val_new(net_sel,j,max_time_idx(j)-bf:max_time_idx(j)+bf);
end
megLI_sub_pt = mean(mLI_sub1,2);


% a = absmax_idx.tsel(time_pts);


if cfg_main.ternary == 1
    cfg = []; cfg.thre = thre;
    cfg.LI = megLI_sub_pt; megLI_sub_pt = do_ternary_classification(cfg);
    Concordance = (length(find([megLI_sub_pt - fmri_LIs_val] == 0))/length(fmri_LIs_val)).*100;
end

disp(['LI concordance (MEG-vs-fMRI): ', num2str(Concordance)])

lgd = {'meg', 'fmri'};
figure, bar([megLI_sub_pt, fmri_LIs_val])
set(gcf, 'Position', [600   500   1000   300]);
set(gca,'color','none');
lgnd = legend(lgd);
% set(lgnd,'color','none','location','east');
set(lgnd,'location','southeast');
L = length(ID);
set(gca,'Xtick', 1:L,'XtickLabel',ID);
set(gca,'FontSize',8,'XTickLabelRotation',90);
axis tight

% - export figs
if savefig == 1
    cfg = [];
    cfg.outdir = outdir;
    cfg.filename = 'MEG_vs_fMRI_concord_subjs';
    cfg.type = 'fig';
    do_export_fig(cfg)
end

%%



end