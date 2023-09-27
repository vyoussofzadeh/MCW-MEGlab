function [megLI_sub_pt, fmri_LIs_val, crr] = do_MEG_fMRI_corr_contrast(cfg_main)


wi = cfg_main.wi;
LI_pt_new = cfg_main.LI_val;
% LI_symb_pt_val_new = cfg_main.LI_symb_val;
net_sel = cfg_main.net_sel;
fmri_LIs_val = cfg_main.fmri_LIs_val;
net_sel_mutiple_label = cfg_main.net_sel_mutiple_label;
ID = cfg_main.ID;
thre = cfg_main.thre;
savefig = cfg_main.savefig;
outdir = cfg_main.outdir;

crr = [];
midx = [];
for i=1:length(wi)
    if length(net_sel) > 1
        
        mLI_sub1 = mean(LI_pt_new(net_sel,:,i));
%         mLI_sub2 = mean(LI_symb_pt_val_new(net_sel,:,i));
    else
        mLI_sub1 = (LI_pt_new(net_sel,:,i));
%         mLI_sub2 = (LI_symb_pt_val_new(net_sel,:,i));
    end
    megLI_sub_pt = (mLI_sub1)';
    
    [~, midx(i)] = max(megLI_sub_pt);
    
    if cfg_main.ternary == 1
        cfg = []; cfg.thre = thre;
        cfg.LI = megLI_sub_pt; megLI_sub_pt = do_ternary_classification(cfg);
    end
    crr(i,:) = corr2(megLI_sub_pt, fmri_LIs_val);
    %     [mx, midx] = max(abs(megLI_sub_pt)); midx_all(i) = midx;
end

% Determine the time point of max LI for each subject
[~, max_time_pts] = max(squeeze(LI_pt_new(net_sel,:,:)), [], 2);
% [~, max_time_pts_symb] = max(squeeze(LI_symb_pt_val_new(net_sel,:,:)), [], 2);

% LI_diff_pt_val = LI_anim_pt_val_new - LI_symb_pt_val_new;
% [val, max_time_pts_diff] = max(squeeze(LI_diff_pt_val(net_sel,:,:)), [], 2);

% difference_LI = squeeze(LI_anim_pt_val_new(net_sel,:,:)) - squeeze(LI_symb_pt_val_new(net_sel,:,:));
% [~, max_time_pts_diff] = max(difference_LI, [], 2);


% figure,
% m1 = squeeze(mean(LI_anim_pt_val_new(net_sel,:,:),1));
% m2 = squeeze(mean(LI_symb_pt_val_new(net_sel,:,:),1));
% m3 = m1 - m2;
% figure, plot(mean(wi'), m3(2,:)), title(ID(2))

figure, plot(mean(wi'),crr), title([net_sel_mutiple_label{net_sel}]);
set(gca,'color','none');
ylabel('LIs corr (MEG vs. fMRI)')
xlabel('Time (sec)')

% - export figs
if savefig == 1
    cfg = [];
    cfg.outdir = outdir;
    cfg.filename = ['corr, ', net_sel_mutiple_label{net_sel}];
    cfg.type = 'fig';
    do_export_fig(cfg)
end

%%
[mx, idx] = max(crr);
bf = cfg_main.bf;

if length(net_sel) > 1
    mLI_sub1 = squeeze(mean(LI_pt_new(net_sel,:,idx-bf:idx+bf)));
%     mLI_sub2 = squeeze(mean(LI_symb_pt_val_new(net_sel,:,idx-bf:idx+bf)));
    megLI_sub_pt = mean((mLI_sub1),2);
else
    mLI_sub1 = mean(LI_pt_new(net_sel,:,idx-bf:idx+bf),3)';
%     mLI_sub2 = mean(LI_symb_pt_val_new(net_sel,:,idx-bf:idx+bf),3)';
    megLI_sub_pt = (mLI_sub1);
end
% megLI_sub_pt = mean((mLI_sub1),2);

%%
% if length(net_sel) > 1
%     mLI_sub1 = arrayfun(@(idx) mean(LI_pt_val_new(net_sel,idx,:)), max_time_pts);
%     mLI_sub2 = arrayfun(@(idx) mean(LI_symb_pt_val_new(net_sel,idx,:)), max_time_pts_symb);
%     
%     %     mLI_sub1 = arrayfun(@(idx) mean(LI_pt_val_new(net_sel,idx,:)), max_time_pts_diff);
%     %     mLI_sub2 = arrayfun(@(idx) mean(LI_symb_pt_val_new(net_sel,idx,:)), max_time_pts_diff);
%     
%     megLI_sub_pt = mean((mLI_sub1 - mLI_sub2),2);
%     megLI_sub_pt = (mLI_sub1);
% else
%     mLI_sub1 = arrayfun(@(idx) LI_pt_val_new(net_sel,idx), max_time_pts);
%     mLI_sub2 = arrayfun(@(idx) LI_symb_pt_val_new(net_sel,idx), max_time_pts_symb);
%     
%     %     mLI_sub1 = arrayfun(@(idx) LI_pt_val_new(net_sel,idx), max_time_pts_diff);
%     %     mLI_sub2 = arrayfun(@(idx) LI_symb_pt_val_new(net_sel,idx), max_time_pts_diff);
%     
%     megLI_sub_pt = (mLI_sub1 - mLI_sub2);
%     megLI_sub_pt = (mLI_sub1);
% end

%%
if cfg_main.ternary == 1
    cfg = []; cfg.thre = thre;
    cfg.LI = megLI_sub_pt; megLI_sub_pt = do_ternary_classification(cfg);
end

if cfg_main.ternary ~= 1
    
    lgd = {'meg', 'fmri'};
    figure, bar([megLI_sub_pt, fmri_LIs_val])
    set(gcf, 'Position', [600   500   1500   300]);
    set(gca,'color','none');
    lgnd = legend(lgd);
    set(lgnd,'color','none');
    xlabel('Subjects')
    L = length(ID);
    set(gca,'Xtick', 1:L,'XtickLabel',ID);
    set(gca,'FontSize',9,'XTickLabelRotation',90);
    
    % - export figs
    if savefig == 1
        cfg = [];
        cfg.outdir = outdir;
        cfg.filename = 'MEG_vs_fMRI_corr_subjs';
        cfg.type = 'fig';
        do_export_fig(cfg)
    end
    
    %     figure, plot(megLI_sub_pt, fmri_LIs_val,'*')
    %     corr2(megLI_sub_pt, fmri_LIs_val)
    %     tbl = table(megLI_sub_pt,fmri_LIs_val);
    %     mdl = fitlm(tbl,'linear');
    %     plot(mdl)
    %     ylabel('LI values, fMRI')
    %     xlabel('LI values, MEG')
    
end

end