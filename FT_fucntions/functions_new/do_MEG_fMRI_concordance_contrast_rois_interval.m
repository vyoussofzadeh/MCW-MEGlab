function [fmri_LIs_val, conc, kappa_time] = do_MEG_fMRI_concordance_contrast_rois_interval(cfg_main)


lang_id = cfg_main.lang_id;
% bf = cfg_main.buffervalue;
net_sel_id = cfg_main.net_sel_id;

wi = cfg_main.wi;
m_wi = mean(wi');

for j = 1:length(lang_id)
    LI_pt_new = cfg_main.LI_val;
    net_sel = net_sel_id(j);
    fmri_LIs_val = cfg_main.fmri_LIs_val;
    net_sel_mutiple_label = cfg_main.net_sel_mutiple_label;
    ID = cfg_main.ID;
    thre = cfg_main.thre;
    savefig = cfg_main.savefig;
    outdir = cfg_main.outdir;
    
    midx = [];
    for i = 1:length(wi)
        if length(net_sel) > 1
            %             mLI_sub1 = mean(LI_pt_new(net_sel,:,i));
            mLI_sub1 = max(LI_pt_new(net_sel,:,i));
        else
            mLI_sub1 = LI_pt_new(net_sel,:,i);
        end
        megLI_sub_pt = (mLI_sub1)';
        
        [~, midx(i)] = max(megLI_sub_pt);
        
        if cfg_main.ternary == 1
            cfg = []; cfg.thre = thre;
            cfg.LI = megLI_sub_pt; mLI_sub_pt_trn = do_ternary_classification2(cfg);
            matches = mLI_sub_pt_trn == fmri_LIs_val;
            numMatches = sum(matches);
            percentageMatch = (numMatches / length(mLI_sub_pt_trn)) * 100;
            conc(j,i,:) = percentageMatch;
        else
            conc(j,i,:) = megLI_sub_pt .* fmri_LIs_val;
        end
        
        % Find the interval with the best concordance
        [mx, idx] = max(conc(j,:));
        interval = wi(idx,:);
        
        % Use the new interval to calculate mLI_sub1
        disp(['selected time for ', cfg_main.lang_id{j}, ':', num2str(interval(1)), '-', num2str(interval(end)), ' ms']);        
        
        MEG_LI = mLI_sub_pt_trn;
        fMRI_LI = fmri_LIs_val;
        % Create a contingency table from the categorical data
        confusion_matrix = crosstab(MEG_LI, fMRI_LI);
        
        % Calculate observed agreement P_o
        total_agreements = sum(diag(confusion_matrix)); % sum of diagonal elements (where the ratings agree)
        total_cases = sum(confusion_matrix, 'all');
        P_o = total_agreements / total_cases;
        
        % Calculate expected agreement P_e
        row_totals = sum(confusion_matrix, 2); % sums of each row
        column_totals = sum(confusion_matrix, 1); % sums of each column
        expected_agreements = sum((row_totals .* column_totals) / total_cases);
        P_e = expected_agreements / total_cases;
        
        % Calculate Cohen's Kappa
        kappa = (P_o - P_e) / (1 - P_e);
        kappa_time(j,i) = kappa;
        
        % Display the result
%         fprintf('Cohen''s Kappa: %f\n', kappa);
    end
end

% figure,
% plot(mean(wi'), kappa_time, 'LineWidth', 3),
% set(gca, 'color', 'none');
% ylabel('LIs concordance (MEG vs. fMRI)')
% xlabel('Time (sec)')
% legend(net_sel_mutiple_label(net_sel_id), 'Location', 'southoutside', 'NumColumns', 5)
% box off
% if isfield(cfg_main, 'title')
%     title(cfg_main.title)
% end

% figure,
% plot(mean(wi'), conc, 'LineWidth', 3),
% set(gca, 'color', 'none');
% ylabel('LIs concordance (MEG vs. fMRI)')
% xlabel('Time (sec)')
% legend(net_sel_mutiple_label(net_sel_id), 'Location', 'southoutside', 'NumColumns', 5)
% box off
% if isfield(cfg_main, 'title')
%     title(cfg_main.title)
% end

plotTimeIntervalData(wi, conc, net_sel_mutiple_label, cfg_main)

if savefig == 1
    cfg = [];
    cfg.outdir = outdir;
    cfg.filename = ['concor, ', net_sel_mutiple_label{net_sel}];
    cfg.type = 'svg';
    do_export_fig(cfg)
end

end