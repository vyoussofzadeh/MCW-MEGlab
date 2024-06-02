function [megLI_sub_pt, fmri_LIs_val, mwi_interval] = do_MEG_fMRI_concordance_contrast(cfg_main)


wi = cfg_main.wi;
LI_pt_val_new = cfg_main.LI_val;
net_sel = cfg_main.net_sel;
fmri_LIs_val = cfg_main.fmri_LIs_val;
net_sel_mutiple_label = cfg_main.net_sel_mutiple_label;
ID = cfg_main.ID;
thre = cfg_main.thre;
savefig = cfg_main.savefig;
outdir = cfg_main.outdir;
bf = cfg_main.buffervalue;

%%
conc = [];
for i=1:length(wi)
    if length(net_sel) > 1
        mLI_sub1 = mean(LI_pt_val_new(net_sel,:,i));
    else
        mLI_sub1 = (LI_pt_val_new(net_sel,:,i));
    end
    megLI_sub_pt = (mLI_sub1)';
    
    if cfg_main.ternary == 1
        cfg = []; cfg.thre = thre;
        cfg.LI = megLI_sub_pt; mLI_sub_pt_trn = do_ternary_classification(cfg);
        matches = mLI_sub_pt_trn == fmri_LIs_val;
        numMatches = sum(matches);
        percentageMatch = (numMatches / length(mLI_sub_pt_trn)) * 100;
        conc(i,:) = percentageMatch; %(mLI_sub_pt_trn .* fmri_LIs_val);
    else
        conc(i,:) = (megLI_sub_pt .* fmri_LIs_val);
    end
end

%%
figure, plot(mean(wi'),mean(conc,2),'LineWidth', 3), title([net_sel_mutiple_label{net_sel}]);
ylabel('LIs conc (MEG , fMRI)')
set(gca,'color','none');
xlabel('Time (sec)')

[mx, idx] = max(mean(conc,2));
% disp(mx)

mwi = mean(wi,2);
mwi_interval = mwi(idx-bf:idx+bf);

% - export figs
if savefig == 1
    cfg = [];
    cfg.outdir = outdir;
    cfg.filename = [net_sel_mutiple_label{net_sel}];
    cfg.type = 'fig';
    do_export_fig(cfg)
end

if length(net_sel) > 1
    mLI_sub1 = squeeze(mean(LI_pt_val_new(net_sel,:,idx-bf:idx+bf)));
    megLI_sub_pt = mean((mLI_sub1),2);
else
    mLI_sub1 = mean(LI_pt_val_new(net_sel,:,idx-bf:idx+bf),3);
    megLI_sub_pt = (mLI_sub1)';
end

if cfg_main.ternary == 1
    cfg = []; cfg.thre = thre;
    cfg.LI = megLI_sub_pt; megLI_sub_pt = do_ternary_classification(cfg);
    Concordance = (length(find([megLI_sub_pt - fmri_LIs_val] == 0))/length(fmri_LIs_val)).*100;
end

disp(['time interval:', num2str(mwi_interval(1)), '-', num2str(mwi_interval(end)), ' sec'])
disp(['LI concordance (MEG-vs-fMRI): ', num2str(Concordance)])

lgd = {'meg', 'fmri'};
figure, bar([megLI_sub_pt, fmri_LIs_val])
set(gcf, 'Position', [600   500   1000   300]);
set(gca,'color','none');
lgnd = legend(lgd);
set(lgnd,'location','southeast');
% set(lgnd,'color','none');
L = length(ID);
set(gca,'Xtick', 1:L,'XtickLabel',ID);
set(gca,'FontSize',8,'XTickLabelRotation',90);

% - export figs
if savefig == 1
    cfg = [];
    cfg.outdir = outdir;
    cfg.filename = 'MEG_vs_fMRI_concord_subjs';
    cfg.type = 'fig';
    do_export_fig(cfg)
end


end