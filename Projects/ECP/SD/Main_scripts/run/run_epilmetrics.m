%% Mutiple linear regression

%- Step: Build table for regression
fMRI_LI  = bestResultsTable.fMRI_LI{3};

% figure, plot(optMEG_LI, fMRI_LI, 'ro')


% Example: pick some columns from final_combined_updt
T = table( ...
    optMEG_LI, ...
    final_combined_updt.Animal_RT, ...
    final_combined_updt.Symbol_RT, ...
    final_combined_updt.Animal_ACC, ...
    final_combined_updt.Symbol_ACC, ...
    final_combined_updt.AEDCount, ...
    final_combined_updt.TLEside, ...
    final_combined_updt.EHQ, ...
    final_combined_updt.NP1WASI_FSIQ, ...
    final_combined_updt.CP_freq, ...
    rsnr_max', ...
    fMRI_LI, ...
    'VariableNames', { ...
        'MEG_LI', ...         % The outcome
        'Animal_RT', ...
        'Symbol_RT', ...
        'Animal_ACC', ...
        'Symbol_ACC', ...
        'AEDCount', ...
        'TLEside', ...
        'EHQ', ...
        'FSIQ' ...            % renamed NP1WASI_FSIQ -> FSIQ for convenience
        'CP_freq', ...
        'rSNR', ...
        'fMRI_LI'
    });

% 3) Convert TLEside to categorical if necessary
if ~iscategorical(T.TLEside)
    T.TLEside = categorical(T.TLEside);
end

% (Optional) Remove rows with any missing values
T(any(ismissing(T),2), :) = [];

% 4) Fit the linear model
lm = fitlm(T, 'MEG_LI ~ Animal_RT + Symbol_RT + Animal_ACC + Symbol_ACC + AEDCount + TLEside + EHQ + FSIQ');
disp(lm);

% 5) Inspect model output:
% - Look at 'lm.Coefficients' to see each predictor's slope, p-value
% - Check 'lm.Rsquared.Ordinary' or 'lm.Rsquared.Adjusted' for model fit

% (Optional) Stepwise model to find best subset of predictors:
predictorList = {'Animal_RT','Symbol_RT','Animal_ACC','Symbol_ACC','AEDCount','TLEside','EHQ','FSIQ'};
lm_stepwise = stepwiselm(T, 'MEG_LI ~ 1','Upper','linear',...
    'PredictorVars', predictorList, 'ResponseVar','MEG_LI');
disp(lm_stepwise);

%% SNR regression
lm_rsnr = fitlm(T, 'rSNR ~ AEDCount + EHQ + CP_freq + FSIQ + ...');
disp(lm_rsnr);

lm_combo = fitlm(T, 'MEG_LI ~ rSNR + AEDCount + EHQ + CP_freq + FSIQ + ...');
disp(lm_combo);

% Pairwise Scatter Matrix
numericVars = {'MEG_LI','fMRI_LI','rSNR', 'Animal_RT', 'Symbol_RT', 'Animal_ACC', 'Symbol_ACC', 'AEDCount','EHQ','CP_freq','FSIQ'};
figure;
plotmatrix(T{:, numericVars});  % transforms table data into numeric

%%
% close all  clc

figure, gscatter(T.EHQ, T.MEG_LI, T.TLEside);      % plots grouped scatter
xlabel('EHQ');                             % x-axis label
ylabel('MEG LI');                          % y-axis label
title('MEG LI vs. EHQ (Grouped by TLEside)');

figure, gscatter(T.Animal_ACC, T.MEG_LI, T.TLEside);      % plots grouped scatter
xlabel('Animal-ACC');                             % x-axis label
ylabel('MEG LI');                          % y-axis label
title('MEG LI vs. Animal_ACC (Grouped by TLEside)');

figure, gscatter(T.EHQ, T.fMRI_LI, T.TLEside);      % plots grouped scatter
xlabel('EHQ');                             % x-axis label
ylabel('fMRI LI');                          % y-axis label
title('fMRI LI vs. EHQ (Grouped by TLEside)');

%%
% 1) Extract numeric columns as a matrix using curly braces
numericT = T{:, {'MEG_LI','rSNR','AEDCount','EHQ','CP_freq','FSIQ'}};

% 2) Ensure TLEside is a valid grouping variable (often categorical)
if ~iscategorical(T.TLEside)
    T.TLEside = categorical(T.TLEside);
end

% 3) Call parallelcoords with the numeric matrix
figure('Color','w');
parallelcoords(numericT, 'Group', T.TLEside);
legend('show');
title('Parallel Coordinates Plot');




