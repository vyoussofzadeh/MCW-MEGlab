function ecpfunc_regression_analysis(bestResultsTable, roiList, final_combined_updt, ~, saveDir)
% ECP_REGRESSION_ANALYSIS
%   Perform multiple regression for MEG vs fMRI LIs, with user-specified covariates.

% 1) For each ROI in "roiList", extract the best method from bestResultsTable
%    then build a single table T with the needed columns from final_combined_updt.

% If you want to handle multiple rois in a loop:
for iROI = 1:length(roiList)
    roiName = roiList{iROI};

    % ----- A) Find row in bestResultsTable for this ROI
    rowIdx = find(strcmp(bestResultsTable.ROI, roiName));
    if isempty(rowIdx)
        warning('ROI %s not found in bestResultsTable.', roiName);
        continue
    end

    % Extract fields
    optMEG_LI      = bestResultsTable.optMEG_LI{rowIdx};
    fMRI_LI        = bestResultsTable.fMRI_LI{rowIdx};
    MEG_fMRI_diff  = optMEG_LI - fMRI_LI;  % simple difference
    MEG_LI_trn     = bestResultsTable.MEG_LI_tern{rowIdx};
    fMRI_LI_trn    = bestResultsTable.fMRI_LI_tern{rowIdx};
    MEG_fMRI_diff_trn = (MEG_LI_trn - fMRI_LI_trn);

    % (Optional) if you computed rSNR or other measures:
    rSNR_L  = bestResultsTable.rSNR_L{rowIdx};
    rSNR_R  = bestResultsTable.rSNR_R{rowIdx};
    % Maybe define rSNR_max, etc.

    % ----- B) Build the big table T using final_combined_updt columns
    %     Check for matching row counts if needed.

    % Some continuous measures from final_combined_updt:
    Animal_RT  = final_combined_updt.Animal_RT; 
    Symbol_RT  = final_combined_updt.Symbol_RT;
    Animal_ACC = final_combined_updt.Animal_ACC;
    Symbol_ACC = final_combined_updt.Symbol_ACC;
    AEDCount   = final_combined_updt.AEDCount;
    EHQ        = final_combined_updt.EHQ;        % continuous
    CP_freq    = final_combined_updt.CP_freq;
    SG_freq    = final_combined_updt.SG_freq;
    LTGTC      = final_combined_updt.LTGTC;      % might be numeric or cat
    TLEside    = final_combined_updt.TLEside;    % categorical (Left/Right/???)
    FSIQ       = final_combined_updt.NP1WASI_FSIQ; % continuous IQ measure
    % Check size:
    nSubj = length(optMEG_LI);
    if length(Animal_RT) ~= nSubj
        warning('Row mismatch for ROI=%s. Skipping.', roiName);
        continue
    end

    % Create a MATLAB table with relevant variables:
    T = table( ...
        optMEG_LI, ...                 % The main outcome
        MEG_fMRI_diff, ...
        MEG_fMRI_diff_trn, ...
        Animal_RT, ...
        Symbol_RT, ...
        Animal_ACC, ...
        Symbol_ACC, ...
        AEDCount, ...
        TLEside, ...
        EHQ, ...
        CP_freq, ...
        SG_freq, ...
        LTGTC, ...
        FSIQ, ...
        fMRI_LI, ...                   % optional predictor
        'VariableNames', { ...
        'optMEG_LI', 'MEG_fMRI_diff','MEG_fMRI_diff_trn',...
        'Animal_RT','Symbol_RT','Animal_ACC','Symbol_ACC',...
        'AEDCount','TLEside','EHQ','CP_freq','SG_freq','LTGTC','FSIQ','fMRI_LI'} );

    % Convert TLEside to categorical if it's not already:
    if ~iscategorical(T.TLEside)
        T.TLEside = categorical(T.TLEside);
    end

    % Optionally remove rows with missing data:
    T = rmmissing(T);  % or T(any(ismissing(T),2),:) = [];

    % ----- C) Check collinearity
    numericVars = {'Animal_RT','Symbol_RT','Animal_ACC','Symbol_ACC','AEDCount','EHQ','CP_freq','SG_freq','FSIQ','fMRI_LI','MEG_fMRI_diff','MEG_fMRI_diff_trn','optMEG_LI'};
    labelShort  = {'AniRT','SymRT','AniACC','SymACC','AED','EHQ','CP','SG','FSIQ','fMRI','MEGfMRIdiff','MEGfMRItrn','MEG_LI'}; 
    % Possibly remove columns that aren't relevant to the correlation check

    plotLowerHalfCorr(T, numericVars, labelShort, sprintf('CorrMatrix_%s', roiName));
    doPlotExport(1, saveDir, sprintf('CorrMatrix_%s', roiName), 'png'); % or 'svg'

    % Possibly examine partial correlation or VIF, etc.

    % ----- D) Fit a regression model
    % Suppose we want to predict 'optMEG_LI' from these variables:
    mdl1 = fitlm(T, 'optMEG_LI ~ Animal_RT + Symbol_RT + Animal_ACC + Symbol_ACC + AEDCount + TLEside + EHQ + CP_freq + SG_freq + FSIQ + fMRI_LI');
    disp(mdl1);

    % Another model:
    mdl2 = fitlm(T, 'MEG_fMRI_diff ~ Animal_RT + Symbol_RT + Animal_ACC + Symbol_ACC + AEDCount + TLEside + EHQ + CP_freq + SG_freq + FSIQ + fMRI_LI');
    disp(mdl2);

    % Summarize or save results
    modelSummary1 = mdl1.Coefficients;
    modelSummary2 = mdl2.Coefficients;
    % Possibly write them to a file or store in a struct

    % ... (any additional analysis or saving)...

end % end for iROI

end % end function
