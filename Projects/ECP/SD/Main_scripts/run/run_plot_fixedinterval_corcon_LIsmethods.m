for i=1:length(LI_method_label)
    
    disp(['correlation analysis: ',LI_method_label{i}])
    
    cfg = []; cfg.wi = wi;
    cfg.ID = sub_MF_pt;
    cfg.thre = MEG_thre;
    cfg.bf = 1;
    cfg.ternary = 0;
    cfg.savefig = 0;
    cfg.outdir = save_dir;
    cfg.net_sel_mutiple_label = net_sel_mutiple_label;
    cfg.net_sel_id = [1,2,6,11];
    cfg.lang_id = {'language_Angular'; 'language_Frontal';'language_Temporal'; 'language_Lateral'};
    cfg.LI_val = LI_pt_val_new.(LI_method_label{i});
    cfg.fmri_LIs_val = fmri_LIs;
    cfg.idx = IB_megfmri;
    cfg.title = LI_method_label{i};
    do_MEG_fMRI_corr_contrast_rois(cfg);
    
    % - export figs
    cfg = []; cfg.outdir = save_dir; filename = ['corr ROIs_', LI_method_label{i}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
    close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

    disp(['Concordance analysis: ',LI_method_label{i}])
    
    disp('---')
    % clc
    cfg = []; cfg.wi = wi;
    cfg.ID = sub_MF_pt;
    cfg.thre = MEG_thre;
    cfg.ternary = 1;
    cfg.savefig = 0;
    cfg.outdir = save_dir;
    cfg.net_sel_mutiple_label = net_sel_mutiple_label;
    cfg.net_sel_id = [1,2,6,11];
    cfg.lang_id = {'language_Angular'; 'language_Frontal';'language_Temporal'; 'language_Lateral'};
    cfg.LI_val = LI_pt_val_new.(LI_method_label{i});
    cfg.fmri_LIs_val = fmri_LIs_trn;
    cfg.idx = IB_megfmri;
    cfg.title = LI_method_label{i};
    cfg.buffervalue = 5;
%     do_MEG_fMRI_concordance_contrast_rois(cfg);
    do_MEG_fMRI_concordance_contrast_rois_interval(cfg);
    
    % - export figs
    cfg = []; cfg.outdir = save_dir; filename = ['concor ROIs_', LI_method_label{i}]; cfg.filename = filename; cfg.type = 'svg'; do_export_fig(cfg)
    close all, combined_path = fullfile(save_dir,[cfg.filename, '.svg']); web(combined_path, '-new');

end

cd(save_dir)















