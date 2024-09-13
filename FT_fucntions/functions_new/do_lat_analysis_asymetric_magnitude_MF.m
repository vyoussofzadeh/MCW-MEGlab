function [LI, LI_max, pow_values] = do_lat_analysis_asymetric_magnitude_MF(cfg_main)

wi = cfg_main.wi;
atlas = cfg_main.atlas;
idx_L = cfg_main.index_L;
idx_R = cfg_main.index_R;
thre = cfg_main.thre;
sinput = cfg_main.sinput;
doavg = cfg_main.doavg;

%- Parcel_based (mean parcels) LI analysis
tmp = load(fullfile(cfg_main.BS_data_dir, sinput)); tmp.ImageGridAmp = tmp.Value;

idx_LR = [idx_L,idx_R];

switch cfg_main.math
    case 'db'
        Fdata = tmp.ImageGridAmp; tidx = tmp.Time < 0; meanBaseline = mean(Fdata(:,tidx),2);
        Fdata = 10 .* log10(abs(bst_bsxfun(@rdivide, Fdata, meanBaseline)));
        tmp.ImageGridAmp = Fdata;
        tmp.ImageGridAmp(tmp.ImageGridAmp < 0) = 0;
    case 'rectif'
        Fdata = tmp.ImageGridAmp;
        Fdata(Fdata < 0) = 0;
        tmp.ImageGridAmp = Fdata;
    case 'rectif_bslnormal'
        Fdata = tmp.ImageGridAmp(idx_LR,:);
        Fdata(Fdata < 0) = 0;
        tidx = tmp.Time < 0; meanBaseline = mean(Fdata(:,tidx),2);
        Fdata = Fdata./meanBaseline;
        tmp.ImageGridAmp(idx_LR,:) = Fdata;       
    case 'rectif_bslnormal_mean'
        
        Fdata_left  = mean(tmp.ImageGridAmp(idx_L,:),1);
        Fdata_right = mean(tmp.ImageGridAmp(idx_R,:),1);
        
        tidx = tmp.Time < 0; 
        
        mbsl_L = mean(Fdata_left(:,tidx),2); Fdata_left = Fdata_left./mbsl_L;
        mbsl_R = mean(Fdata_right(:,tidx),2); Fdata_right = Fdata_right./mbsl_R;
        
        for j = 1:length(idx_L), tmp.ImageGridAmp(idx_L(j),:) = Fdata_left; end
        for j = 1:length(idx_R), tmp.ImageGridAmp(idx_R(j),:) = Fdata_right; end
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
    cfg.parcellaion = cfg_main.parcellaion;
    cfg.globalmax = max(max(tmp.ImageGridAmp));
    switch cfg_main.math
        case 'rectif_bslnormal_mean'
            cfg.applymean = 1;
        otherwise
            cfg.applymean = 0;    
    end

    [LI_clin, pow] = do_LI_magnitude(cfg);
    LI(j) = LI_clin;
    
    % Store each pow struct in the pow_values array
    pow_values.left = [pow_values.left; pow.left];
    pow_values.right = [pow_values.right; pow.right];
end

% pow_values.left_org = mean(tmp.ImageGridAmp(idx_L,:),1)';
% pow_values.right_org = mean(tmp.ImageGridAmp(idx_R,:),1)';

if cfg_main.fplot ==1
    figure;
    yyaxis left; % Left y-axis for LI
    plot(LI,'LineWidth',1.5);
    ylabel('Lateralization Index (LI)');
    
    hold on;
    
    yyaxis right; % Right y-axis for Power
    plot(mean(pow_values.left,2),'LineWidth',1.5); % Mean power for left activities
    plot(mean(pow_values.right,2),'LineWidth',1.5); % Mean power for right activities
    ylabel('Power');
    
    legend({'LI', 'Left SMag', 'Right SMag'}, 'Location', 'best');
    
    % Set x-axis ticks and labels
    val = round(mean(wi(:,1),2),2);
    set(gca, 'Xtick', 1:2:length(wi), 'XtickLabel', val(1:2:end));
    set(gca, 'FontSize', 8, 'XTickLabelRotation', 90);
    set(gcf, 'Position', [1000, 400, 1000, 300]);
    
    xlabel('Mean Temporal Windows (sec)');
    title(['SMag']);
    set(gca, 'color', 'none'); % Transparent background
    
    
    figure, plot(tmp.Time, tmp.ImageGridAmp(idx_L,:),'LineWidth',1.5); title('Left H');% Mean power for left activities
    xlim([tmp.Time(1), tmp.Time(end)])
    figure, plot(tmp.Time, tmp.ImageGridAmp(idx_R,:),'LineWidth',1.5); title('Right H');% Mean power for left activities
    xlim([tmp.Time(1), tmp.Time(end)])

end

[~, idx_mx] = max(LI); LI_max = wi(idx_mx,:);
