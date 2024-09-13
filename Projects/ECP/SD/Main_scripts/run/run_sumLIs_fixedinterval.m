% Initialize a table to store results
resultsTable = table([], [], 'VariableNames', {'Method', 'Metrics'});

for i=1:length(LI_method_label)
    
    disp(['Processing: ', LI_method_label{i}]);
    
    % Configuration for both analyses
    cfg = [];
    cfg.wi = wi;
    cfg.ID = sub_MF_pt;
    cfg.savefig = 0; % Assuming this controls figure saving inside the functions
    cfg.outdir = save_dir;
    cfg.net_sel_mutiple_label = net_sel_mutiple_label;
    cfg.net_sel_id = [1,2,6,11];
    cfg.lang_id = {'language_Angular'; 'language_Frontal';'language_Temporal'; 'language_Lateral'};
    cfg.LI_val = LI_pt_val_new.(LI_method_label{i});
    cfg.fmri_LIs_val = fmri_LIs;
    cfg.idx = IB_megfmri;
    cfg.title = LI_method_label{i};
    
    % Correlation Analysis
    cfg.thre = MEG_thre;
    cfg.ternary = 0;
    [megLI_sub, fmri_LIs, correlationMetrics] = do_MEG_fMRI_corr_contrast_rois2(cfg);
    
    % Concordance Analysis
    cfg.thre = MEG_thre;
    cfg.ternary = 1;
    cfg.buffervalue = 2;
    cfg.fmri_LIs_val = fmri_LIs_trn;
%     [~, ~, concordanceMetrics] = do_MEG_fMRI_concordance_contrast_rois(cfg);
    [~, concordanceMetrics] = do_MEG_fMRI_concordance_contrast_rois_interval(cfg);
    
    % Compile results
    metrics = struct('Correlation', correlationMetrics, 'Concordance', concordanceMetrics);
    resultsTable = [resultsTable; {LI_method_label{i}, metrics}];
end
close all,
