function [megLI_sub_pt, fmri_LIs_val, crr_all] = do_MEG_fMRI_corr_contrast_rois2(cfg_main)


% lang_id = {'language_Angular'; 'language_Frontal'; 'language_Occipital'; 'language_Other'; 'language_PCingPrecun'; 'language_Temporal'; 'language_Lateral'};

lang_id = cfg_main.lang_id;

net_sel_id = cfg_main.net_sel_id;
% net_sel_id = [1,2,3,4,5,6,11];


crr_all = [];
anov_F_all = [];

for j=1:length(lang_id)
    
    wi = cfg_main.wi;
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
    anov_F =[];
    
    for i=1:length(wi)
        if length(net_sel) > 1
            mLI_sub1 = nanmean(LI_pt_new(net_sel,:,i));
        else
            mLI_sub1 = (LI_pt_new(net_sel,:,i));
        end
        megLI_sub_pt = (mLI_sub1)';
        
        [~, midx(i)] = max(megLI_sub_pt);
        
        if cfg_main.ternary == 1
            cfg = []; cfg.thre = thre;
            cfg.LI = megLI_sub_pt; megLI_sub_pt = do_ternary_classification2(cfg);
        end
        tmp = fmri_LIs_val.val.(lang_id{j});
        fMRILI_sub_pt = tmp(cfg_main.idx);
        
        %% anova test (optional)
        % Fit the repeated measures model for Source Magnitude
        t_SM = table(megLI_sub_pt, fMRILI_sub_pt, 'VariableNames', {'MEG_SM', 'fMRI_SM'});
        % Fit the repeated measures model for Source Magnitude
        rm_SM = fitrm(t_SM, 'MEG_SM-fMRI_SM~1');
        ranova_SM = ranova(rm_SM);        
        anov_F(i) = ranova_SM.F(1);

        %% Correlation (non-parametric)
        % Calculate Spearman's rank correlation coefficient
        [rho, pval] = corr(megLI_sub_pt, fMRILI_sub_pt, 'Type', 'Spearman');
        
        % Display the correlation coefficient and p-value
%         fprintf('Spearman''s rho: %f\n', rho);
%         fprintf('P-value: %f\n', pval);
        
        crr(i,:) = rho;

        %% Correlation (parametric)
%         [groupCorrelation, pval] = corr(megLI_sub_pt, fMRILI_sub_pt, 'Rows', 'complete');
%         crr(i,:) = groupCorrelation;
        %         crr(i,:) = corr2(megLI_sub_pt, fMRILI_sub_pt);
    end
    crr_all(j,:) = crr;
    anov_F_all(j,:) = anov_F;
end

% figure,
% plot(wi(:,1)',anov_F_all,'LineWidth', 3),
% % title([net_sel_mutiple_label{net_sel}]);
% set(gca,'color','none');
% ylabel('Anova F (MEG vs. fMRI)')
% xlabel('Time (sec)')
% legend(net_sel_mutiple_label(net_sel_id),'Location','southoutside', 'NumColumns', 5)
% box off
% if isfield(cfg_main, 'title')
%     title(cfg_main.title)
% end

figure,
plot(wi(:,1)',crr_all,'LineWidth', 3),
% title([net_sel_mutiple_label{net_sel}]);
set(gca,'color','none');
ylabel('LIs corr (MEG vs. fMRI)')
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
    cfg.filename = ['corr, ', net_sel_mutiple_label{net_sel}];
    cfg.type = 'svg';
    do_export_fig(cfg)
    combined_path = fullfile(outdir,[cfg.filename, '.svg']); web(combined_path, '-new');
end


end