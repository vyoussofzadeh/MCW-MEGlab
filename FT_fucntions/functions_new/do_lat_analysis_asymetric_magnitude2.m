function [LI, LI_max, pow_values] = do_lat_analysis_asymetric_magnitude2(cfg_main)

wi = cfg_main.wi;
atlas = cfg_main.atlas;
idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;
doavg = cfg_main.doavg;


thre = cfg_main.thre;
sinput = cfg_main.sinput;

% tmp = load(fullfile(cfg_main.BS_data_dir, sinput));

tmp_1 = load(fullfile(cfg_main.BS_data_dir, sinput{1}));
tmp_2 = load(fullfile(cfg_main.BS_data_dir, sinput{2}));

% disp({tmp_1.Comment; tmp_2.Comment})

tmp = tmp_1;

switch cfg_main.math
    case 'db'
        tmp.ImageGridAmp = 10 * log10(tmp_1.ImageGridAmp) - 10 * log10(tmp_2.ImageGridAmp); tmp.ImageGridAmp(tmp.ImageGridAmp < 0) = 0;
        
    case 'db_baseline'
        Fdata1 = tmp_1.ImageGridAmp; tidx = tmp_1.Time < 0; meanBaseline = mean(Fdata1(:,tidx),2);
        Fdata1 = 10 .* log10(abs(bst_bsxfun(@rdivide, Fdata1, meanBaseline)));
        
        Fdata2 = tmp_2.ImageGridAmp; tidx = tmp_2.Time < 0; meanBaseline = mean(Fdata2(:,tidx),2);
        Fdata2 = 10 .* log10(abs(bst_bsxfun(@rdivide, Fdata2, meanBaseline)));
        
        tmp.ImageGridAmp = Fdata1 - Fdata2;
        tmp.ImageGridAmp(tmp.ImageGridAmp < 0) = 0;
        
    case 'diff'
        tmp.ImageGridAmp = tmp_1.ImageGridAmp - tmp_2.ImageGridAmp;
        
    case 'rec_diff'
        tmp.ImageGridAmp = tmp_1.ImageGridAmp - tmp_2.ImageGridAmp;
        tmp.ImageGridAmp(tmp.ImageGridAmp < 0) = 0;
        
    case 'rec_diff_mbsl'
        tmp.ImageGridAmp(idx_LR_updt,:) = tmp_1.ImageGridAmp(idx_LR_updt,:) - tmp_2.ImageGridAmp(idx_LR_updt,:);
        Fdata = tmp.ImageGridAmp(idx_LR_updt,:); Fdata(Fdata < 0) = 0;
        tidx = tmp.Time < 0;
        meanBaseline = mean(Fdata(:,tidx),2);
        Fdata = Fdata./meanBaseline; tmp.ImageGridAmp(idx_LR_updt,:) = Fdata;
        %         figure, plot(tmp.Time, Fdata');

    case 'rec_diff_zbsl'
        tmp.ImageGridAmp(idx_LR_updt,:) = tmp_1.ImageGridAmp(idx_LR_updt,:) - tmp_2.ImageGridAmp(idx_LR_updt,:);
        Fdata = tmp.ImageGridAmp(idx_LR_updt,:); Fdata(Fdata < 0) = 0;
        tidx = tmp.Time < 0;
        meanBaseline = mean(Fdata(:,tidx),2);
        stdBaseline = std(Fdata(:,tidx)')';
        Fdata = (Fdata - meanBaseline)./stdBaseline; tmp.ImageGridAmp(idx_LR_updt,:) = Fdata;       

    case 'bsl_diff_rec'
        Fdata = tmp_1.ImageGridAmp(idx_LR_updt,:);
        tidx = tmp.Time < 0; meanBaseline = mean(Fdata(:,tidx),2); Fdata = Fdata./meanBaseline; % bsl norm
        tmp_1.ImageGridAmp(idx_LR_updt,:) = Fdata;
        
        Fdata = tmp_2.ImageGridAmp(idx_LR_updt,:);
        tidx = tmp.Time < 0; meanBaseline = mean(Fdata(:,tidx),2); Fdata = Fdata./meanBaseline; % bsl norm
        tmp_2.ImageGridAmp(idx_LR_updt,:) = Fdata;
        
        tmp.ImageGridAmp(idx_LR_updt,:) = tmp_1.ImageGridAmp(idx_LR_updt,:) - tmp_2.ImageGridAmp(idx_LR_updt,:);
        Fdata = tmp.ImageGridAmp(idx_LR_updt,:);
        Fdata(Fdata < 0) = 0;
        tmp.ImageGridAmp(idx_LR_updt,:) = Fdata;
end

LI = []; 
% Initialize an array to store pow values for all intervals
pow_values = struct('left', [], 'right', []); % Create a struct to hold all pow values

for j=1:size(wi,1)
    
    timind1 = nearest(tmp.Time, wi(j,1)); timind2 = nearest(tmp.Time, wi(j,2));
    
    cfg = [];
    cfg.thre = thre;
    cfg.atlas = atlas;
    if doavg == 1
        cfg.d_in = mean(tmp.ImageGridAmp(:,timind1:timind2),2);
    else
        cfg.d_in = tmp.ImageGridAmp(:,timind1:timind2);
    end

    cfg.idx_L = idx_L;
    cfg.idx_R = idx_R;
    cfg.Threshtype = cfg_main.Threshtype;
    cfg.thre = thre;
    cfg.globalmax = max(max(tmp.ImageGridAmp));
    [LI_clin, pow] = do_LI_magnitude(cfg);
    LI(j) = LI_clin;
    
    % Store each pow struct in the pow_values array
    pow_values.left = [pow_values.left; pow.left];
    pow_values.right = [pow_values.right; pow.right];
end

if cfg_main.fplot ==1
    figure,plot(LI),
    val = round(mean(wi(:,1),2),2);
    set(gca,'Xtick', 1:2:length(wi),'XtickLabel',val(1:2:end));
    set(gca,'FontSize',8,'XTickLabelRotation',90);
    set(gcf, 'Position', [1000   400   1000   300]);
    xlabel('temporal windows (sec)')
    ylabel('LI')
    set(gca,'color','none');
end

[~, idx_mx] = max(LI); LI_max = wi(idx_mx,:);

end
