function [megLI_sub_pt, fmri_LIs_val, crr] = do_MEG_fMRI_concordance_contrast_rois(cfg_main)


% lang_id = {'language_Angular'; 'language_Frontal'; 'language_Occipital'; 'language_Other'; 'language_PCingPrecun'; 'language_Temporal'; 'language_Lateral'};

lang_id = cfg_main.lang_id;
bf = cfg_main.buffervalue;
net_sel_id = cfg_main.net_sel_id;
LI_pt_val_new1 = cfg_main.LI_val;
% net_sel_id = [1,2,3,4,5,6,11];

wi = cfg_main.wi;
m_wi = mean(wi');


conc = [];
for j=1:length(lang_id)
    
    
    LI_pt_new = cfg_main.LI_val;
    net_sel = net_sel_id(j);
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
        else
            mLI_sub1 = (LI_pt_new(net_sel,:,i));
        end
        megLI_sub_pt = (mLI_sub1)';
        
        [~, midx(i)] = max(megLI_sub_pt);
        
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
    
    [mx, idx] = max(conc(j,:));
    interval = idx-bf:idx+bf;
    
    mLI_sub1 = mean(LI_pt_new(net_sel,:,interval),3);
    megLI_sub_pt = (mLI_sub1)';
    disp(['selected time for ', cfg_main.lang_id{j}, ':', num2str(m_wi(interval(1))),'-',num2str(m_wi(interval(end))),' ms'])

end

figure,
plot(mean(wi'),conc,'LineWidth', 3),
% title([net_sel_mutiple_label{net_sel}]);
set(gca,'color','none');
ylabel('LIs concordance (MEG vs. fMRI)')
xlabel('Time (sec)')
legend(net_sel_mutiple_label(net_sel_id),'Location','southoutside', 'NumColumns', 5)
box off
if isfield(cfg_main, 'title')
    title(cfg_main.title)
end

% - export figs
if savefig == 1
    cfg = [];
    cfg.outdir = outdir;
    cfg.filename = ['concor, ', net_sel_mutiple_label{net_sel}];
    cfg.type = 'fig';
    do_export_fig(cfg)
end

end