

%% === Firth logistic regression: Discordant (1) vs Non-discordant (0) ===
% Outcome vector
y = zeros(height(T),1);
y(discordSubs) = 1;

% Candidate predictors (present in T). Mix of continuous & categorical.
candCont = {'Animal_RT','Symbol_RT','Animal_ACC','Symbol_ACC', ...
            'AEDCount','EHQ','CP_freq','SG_freq','NP1WASI_FSIQ','rSNR', ...
            'nSNR_Beta_tSSS_vs_Raw','nSNR_Beta_MEGnet_vs_tSSS', ...
            'nSNR_Broad_tSSS_vs_Raw','nSNR_Broad_MEGnet_vs_tSSS'};
candCat  = {'TLEside','LTGTC'};  % add others if available in T

candCont = candCont(ismember(candCont, T.Properties.VariableNames));
candCat  = candCat(ismember(candCat,  T.Properties.VariableNames));

% --- Univariate screens (logit) & BH-FDR ranking ---
varNames = [candCont, candCat];
p_uni    = nan(numel(varNames),1);

for k = 1:numel(varNames)
    vn = varNames{k};
    x  = T.(vn);

    tbl = table(y);
    if iscell(x) || isstring(x) || iscategorical(x)
        x = categorical(x);
        tbl.x = x;
        mdl   = fitglm(tbl, 'y ~ x', 'Distribution','binomial', 'CategoricalVars','x');
    else
        x = double(x);
        % z-scale (omit NaNs)
        mu = mean(x,'omitnan'); sd = std(x,'omitnan');
        if sd==0 || isnan(sd), sd = 1; end
        xz = (x - mu) ./ sd;
        tbl.xz = xz;
        mdl    = fitglm(tbl, 'y ~ xz', 'Distribution','binomial');
    end

    % overall effect of this predictor (handles multi-level categorical)
    p_uni(k) = coefTest(mdl);   % ?² test that all non-intercept betas = 0
end

% BH-FDR (BenjaminiHochberg)
[~, sortIdx] = sort(p_uni);                          % ascending p
m = sum(~isnan(p_uni));
q = nan(size(p_uni));
ranks = zeros(size(p_uni));
ranks(sortIdx) = 1:numel(sortIdx);
q(sortIdx) = p_uni(sortIdx) .* m ./ ranks(sortIdx);
% make monotone
for i = numel(sortIdx)-1:-1:1
    q(sortIdx(i)) = min(q(sortIdx(i)), q(sortIdx(i+1)));
end

% Pick top 3 by FDR (break ties by raw p)
valid = ~isnan(q);
[~, ord] = sortrows([q(valid), p_uni(valid)], [1 2]);
topIdx_all = find(valid);
topIdx = topIdx_all(ord(1:min(3, numel(ord))));
topVars = varNames(topIdx);

fprintf('\nTop covariates by FDR (q):\n');
for ii = 1:numel(topVars)
    fprintf('  %s  (p=%.3g, q=%.3g)\n', topVars{ii}, p_uni(topIdx(ii)), q(topIdx(ii)));
end

% --- Build design matrix: z-scale continuous, dummy-code categoricals ---
rowMask = ~isnan(y);
X = []; namesX = {};

% precompute masks for missingness per selected var
for k = 1:numel(topVars)
    vn = topVars{k};
    x  = T.(vn);
    if iscell(x) || isstring(x) || iscategorical(x)
        rowMask = rowMask & ~ismissing(x);
    else
        rowMask = rowMask & ~isnan(x);
    end
end

% construct X (no intercept here; add later consistently for both solvers)
for k = 1:numel(topVars)
    vn = topVars{k};
    x  = T.(vn);

    if iscell(x) || isstring(x) || iscategorical(x)
        x = categorical(x);
        D = dummyvar(x(rowMask));          % k levels -> k columns
        if size(D,2) > 1                   % drop first as reference
            D = D(:,2:end);
            cats = categories(x);
            for j = 2:numel(cats)
                namesX{end+1} = sprintf('%s_%s', vn, cats{j}); %#ok<AGROW>
            end
        else
            namesX{end+1} = sprintf('%s_%s', vn, char(categories(x)));
        end
        X = [X, D]; %#ok<AGROW>
    else
        x = double(x);
        mu = mean(x(rowMask),'omitnan'); sd = std(x(rowMask),'omitnan'); if sd==0, sd=1; end
        xz = (x(rowMask) - mu) ./ sd;      % z-scale
        X  = [X, xz];                      %#ok<AGROW>
        namesX{end+1} = ['z_' vn];         %#ok<AGROW>
    end
end

%%
ySel = y(rowMask);

% --- Fit Firth logistic if available; else standard logistic as fallback ---
useFirth = (exist('firthlogit','file')==2);
if useFirth
    % Many firthlogit functions expect X WITHOUT intercept; add internally if needed.
    try
        [beta, se, stats] = firthlogit([ones(size(X,1),1) X], ySel);   % (Intercept first)
        b     = beta(:);
        bSE   = se(:);
        pvals = stats.p(:);
    catch
        % Alternative signature: firthlogit(X,y) assumes you add intercept yourself
        [beta, se, stats] = firthlogit(X, ySel);
        % prepend intercept by re-fitting with an intercept column via augmentation
        % (skip if your firthlogit already includes it)
        b     = beta(:);
        bSE   = se(:);
        pvals = stats.p(:);
    end
    % If a method-specific CI not provided, use Wald as a reasonable approximation
    z = 1.96;
    OR      = exp(b);
    OR_CI_L = exp(b - z*bSE);
    OR_CI_U = exp(b + z*bSE);

    % Build result table (drop intercept if present)
    if numel(b) == numel(namesX) + 1
        varOut = [{'Intercept'}, namesX];
    else
        varOut = namesX;
    end
    res = table(varOut(:), b(:), bSE(:), pvals(:), OR(:), OR_CI_L(:), OR_CI_U(:), ...
        'VariableNames', {'Term','Beta','SE','p','OR','OR_L','OR_U'});
else
    warning('Firth logistic not found on path. Falling back to standard logistic (fitglm).');
    tblGLM = array2table(X, 'VariableNames', namesX);
    tblGLM.y = ySel;
    mdl = fitglm(tblGLM, 'y ~ 1 + ...', 'Distribution','binomial');  % includes intercept
    co = mdl.Coefficients;
%     res = table(co.Properties.RowNames, co.Estimate, co.SE, co.pValue, ...
%                 exp(co.Estimate), exp(co.Lower), exp(co.Upper), ...
%                 'VariableNames', {'Term','Beta','SE','p','OR','OR_L','OR_U'});
end
% 
% disp('--- Multivariable (Firth) logistic results ---');
% disp(res);
% 
% % Optional: save results
% try
%     save(fullfile(save_dir, sprintf('FirthLogit_%s.mat', roi)), 'res', 'topVars', 'q', 'p_uni', 'rowMask', 'namesX');
% catch, end


%%
% After: mdl = fitglm(tblGLM, formulaStr, 'Distribution','binomial');

co  = mdl.Coefficients;
dfe = mdl.DFE;  % degrees of freedom

% Replace NaN p-values using t (if dfe>0) else normal (z) approximation
p_from_t  = 2*tcdf(-abs(co.tStat), max(dfe,1));   % safe even if dfe==0
p_from_z  = 2*normcdf(-abs(co.tStat));            % large-sample approx
p_fix     = p_from_t;
if ~isfinite(dfe) || dfe<=0
    p_fix = p_from_z;
end
nanMask = isnan(co.pValue);
co.pValue(nanMask) = p_fix(nanMask);

% 95% CI for betas (use coefCI if available; else Wald)
try
    CI = coefCI(mdl);   % [L U] for each term
catch
    z = 1.96;
    CI = [co.Estimate - z*co.SE, co.Estimate + z*co.SE];
end

% Odds ratios and CIs
OR    = exp(co.Estimate);
OR_L  = exp(CI(:,1));
OR_U  = exp(CI(:,2));

% Pack results
res = table( ...
    string(co.Properties.RowNames), co.Estimate, co.SE, co.tStat, co.pValue, ...
    OR, OR_L, OR_U, ...
    'VariableNames', {'Term','Beta','SE','Stat','p','OR','OR_L','OR_U'});

disp(res)

% Optional: sanity checks
fprintf('GLM: N=%d, p=%d, DFE=%d\n', mdl.NumObservations, mdl.NumEstimatedCoefficients, mdl.DFE);

%%
% AUC
scores = predict(mdl, tblGLM);
[~,~,~,auc] = perfcurve(ySel, scores, 1);
fprintf('AUC = %.3f\n', auc);

% McFadden's R^2
mdl0 = fitglm(tblGLM, 'y ~ 1', 'Distribution','binomial');
mcf = 1 - mdl.LogLikelihood / mdl0.LogLikelihood;
fprintf('McFadden R^2 = %.3f\n', mcf);

%%
Xf = tblGLM{:, namesX};   % predictors (no intercept)
yf = ySel;                 % 0/1 outcome (discordant)
res_firth = firthlogit_simple(Xf, yf, ["Intercept", string(namesX)]);
disp(res_firth)
