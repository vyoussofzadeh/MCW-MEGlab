function [megLI_sub_pt_trn_All, fmri_LIs_val] = do_MEG_fMRI_concordance_contrast_approches(cfg_main)


wi = cfg_main.wi;
LI_pt_val_new1 = cfg_main.LI_val;
net_sel = cfg_main.net_sel;
fmri_LIs_val = cfg_main.fmri_LIs_val;
net_sel_mutiple_label = cfg_main.net_sel_mutiple_label;
ID = cfg_main.ID;
thre = cfg_main.thre;
savefig = cfg_main.savefig;
outdir = cfg_main.outdir;
bf = cfg_main.buffervalue;
LI_method_label = cfg_main.LI_method_label;

%%
m_wi = mean(wi');

conc = [];
for j=1:size(LI_method_label,2)
    LI_pt_val_new = LI_pt_val_new1.(LI_method_label{j});
    
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
            conc(j,i,:) = percentageMatch; %(mLI_sub_pt_trn .* fmri_LIs_val);
        else
            conc(j,i,:) = (megLI_sub_pt .* fmri_LIs_val);
        end
    end
%     mLI_sub1 = [];
    [mx, idx] = max(conc(j,:)); 
    interval = idx-bf:idx+bf;
    if length(net_sel) > 1
        mLI_sub1 = squeeze(nanmean(LI_pt_val_new(net_sel,:,interval)));
        megLI_sub_pt = nanmean((mLI_sub1),2);
    else
        mLI_sub1 = mean(LI_pt_val_new(net_sel,:,interval),3);
        megLI_sub_pt = (mLI_sub1)';
    end
    disp(['selected time for ', LI_method_label{j}, ':', num2str(m_wi(interval(1))),'-',num2str(m_wi(interval(end))),' ms'])
    
    if cfg_main.ternary == 1
        cfg = []; cfg.thre = thre;
        cfg.LI = megLI_sub_pt; megLI_sub_pt_trn = do_ternary_classification(cfg);
%         (length(find([megLI_sub_pt_trn - fmri_LIs_val] == 0))/length(fmri_LIs_val)).*100
        megLI_sub_pt_trn_All{j} = megLI_sub_pt_trn;
    end
end

%%
figure, plot(mean(wi'),conc,'LineWidth', 3), title([net_sel_mutiple_label{net_sel}]);
ylabel('LIs conc (MEG , fMRI)')
set(gca,'color','none');
xlabel('Time (sec)')
box off
legend(LI_method_label)

[mx, ~] = max(conc');
disp(mx)

% - export figs
if savefig == 1
    cfg = [];
    cfg.outdir = outdir;
    cfg.filename = [net_sel_mutiple_label{net_sel}];
    cfg.type = 'fig';
    do_export_fig(cfg)
end


end