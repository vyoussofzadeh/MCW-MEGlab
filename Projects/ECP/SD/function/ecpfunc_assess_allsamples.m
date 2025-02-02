function ecpfunc_assess_allsamples(cfg_main)

%% Setup / Initialization
roi_sel = cfg_main.roi_sel;
bestResultsTable = cfg_main.bestResultsTable;

discordSubs  = bestResultsTable.Best_Discord_Subs{roi_sel};
roi          = bestResultsTable.ROI{roi_sel};
method       = bestResultsTable.Best_LI_Method{roi_sel};
MEG_LI       = bestResultsTable.MEG_LI{roi_sel};
fMRI_LI      = bestResultsTable.fMRI_LI{roi_sel};
rSNR_L       = bestResultsTable.rSNR_L{roi_sel};
rSNR_R       = bestResultsTable.rSNR_R{roi_sel};
optMEG_LI    = bestResultsTable.optMEG_LI{roi_sel};


wi = cfg_main.wi;
bounds = cfg_main.bounds;
plot_option = cfg_main.plot_option;
save_dir = cfg_main.save_dir;
final_combined_updt = cfg_main.final_combined_updt;
T1_epil_measures = cfg_main.T1_epil_measures;
Pt_ID = final_combined_updt.SubjectID;

Animal_ACC = final_combined_updt.Animal_ACC;
Symbol_ACC = final_combined_updt.Symbol_ACC;

Animal_RT = final_combined_updt.Animal_RT;
Symbol_RT = final_combined_updt.Symbol_RT;

%%
subIDs    = 1:72;  % or however you store subject labels

figure('Color','w','Name','MEG vs fMRI LI','Position',[200,200,700,500]);
scatter(optMEG_LI, fMRI_LI, 50, 'ro', 'filled');  

hold on; grid on;
xlabel('Optimal MEG LI');
ylabel('fMRI LI');
title('MEG LI vs. fMRI LI');

% Add Subject Labels
for s = 1:length(optMEG_LI)
    text(optMEG_LI(s), fMRI_LI(s), num2str(subIDs(s)), ...
        'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
        'FontSize',8, 'Color','k');
end

xLimits = xlim;
plot(xLimits, xLimits, 'k--', 'LineWidth',1);

% (Optional) Compute and Display Correlation / Regression
[r, p] = corr(optMEG_LI, fMRI_LI, 'Rows','complete');
% Display in the plot:
corrStr = sprintf('r = %.2f (p=%.3g)', r, p);
text(xLimits(1) + 0.05*range(xLimits), ...
     xLimits(2) - 0.1*range(xLimits), ...
     corrStr, 'FontSize',10, 'Color','b');

% Or do a quick linear fit:
pCoeffs = polyfit(optMEG_LI, fMRI_LI, 1);  % slope/intercept
xVals = linspace(min(optMEG_LI), max(optMEG_LI), 100);
yFit  = polyval(pCoeffs, xVals);
plot(xVals, yFit, 'b-', 'LineWidth',1.5);

legend({'Data Points','y=x','Corr line'}, 'Location','best');

% Styling Adjustments
% - Adjust axis tightness, fonts, etc.
axis tight;
set(gca,'FontName','Helvetica','FontSize',10);

%% Scatter Plots for Each Covariate vs. MEG_LI
% Suppose T is the table you created with columns:
%   T.MEG_LI     (the measure you want to compare)
%   T.AEDCount   (epilepsy measure)
%   T.EHQ        (epilepsy measure)
%   T.CP_freq    (epilepsy measure)
%   ... etc.

% 1) List your epilepsy covariates in an array
% 'TLEside'
epilepsyCovars = {'AEDCount','EHQ','CP_freq','FSIQ', 'EHQ', 'FSIQ'};  % Example measures

% 2) Optionally define subject IDs for labeling
subIDs = 1:height(T);  % or a cell array of actual subject strings if you have them

% 3) Loop through each covariate and make a scatter plot
figure('Name','Epilepsy Covariates vs MEG LI','Color','w','Position',[100,100,1200,600]);

numCovars = length(epilepsyCovars);
for iCov = 1:numCovars
    subplot(2, ceil(numCovars/2), iCov);  % arrange subplots
    xData = T.(epilepsyCovars{iCov});
    yData = T.MEG_LI;
    yData = T.rsnr_max;
    yData = T.fMRI_LI
    
    scatter(xData, yData, 40, 'b', 'filled');
    hold on; grid on;
    xlabel(epilepsyCovars{iCov});
    ylabel('MEG-LI');
    title(['MEG-LI vs ', epilepsyCovars{iCov}]);
    
    % Add subject labels if desired
%     for s = 1:height(T)
%         text(xData(s), yData(s), num2str(subIDs(s)), ...
%             'VerticalAlignment','bottom', 'HorizontalAlignment','left', ...
%             'FontSize',7, 'Color','k');
%     end
    
    % Add correlation line or stats
    [r, p] = corr(xData, yData, 'Rows','complete');
    corrStr = sprintf('r=%.2f, p=%.3g', r, p);
    xLimits = xlim;
    yLimits = ylim;
    text(xLimits(1)+0.05*range(xLimits), yLimits(2)-0.1*range(yLimits), ...
        corrStr, 'Color','r', 'FontSize',8);

    % Optional: regression line
    pCoeffs = polyfit(xData, yData, 1);
    xVals   = linspace(min(xData), max(xData), 100);
    yFit    = polyval(pCoeffs, xVals);
    plot(xVals, yFit, 'r-', 'LineWidth',1.5);
    
    hold off;
end


end

