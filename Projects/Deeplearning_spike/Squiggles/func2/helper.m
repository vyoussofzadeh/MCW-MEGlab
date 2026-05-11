% Use TRAIN set to avoid leakage (after normalization + fixed T)
fs = 1000;              % put your sampling rate here
r  = 306;                % channel PCs to keep
Ssubj = pca_timeseries_per_subject(XTrain2, YTrain, groups(trainMask), r, fs);

% If you want only a few subjects:
Ssubj = pca_timeseries_per_subject(XTrain2, YTrain, groups(trainMask), r, fs, {'alby'});


S = get_subject_epochs(trials, 'alby');
Xspk  = S.XSpike;      % {nSpike×1} cell, each C×T
Xnspk = S.XNoSpike;    % {nNoSpike×1} cell
disp([S.subj, '  spikes=', num2str(S.nSpike), '  nospike=', num2str(S.nNoSpike)])

D = [];
for i=1:length(Xspk)
    D(:,:,i) = Xspk{i, 1}(:,1:67);
end

figure, plot(mean(D,3)')

cfg = [];
cfg.blocksize = 10;
cfg.viewmode = 'vertical'; %butterfly';
cfg.continuous = 'yes';
cfg.axisfontsize = 7;
cfg.fontsize = 7;
cfg.channel = 'MEG*';
%     cfg.preproc.demean = 'yes';
cfg.position = [300   900   500   1500];
ft_databrowser(cfg, Xspk{1, 1});

%%
separability_pcaLDA(XTrain2, YTrain, groups(trainMask), 16, 0);
separability_CSP(XTrain2, YTrain, groups(trainMask), 3, 0);

evalFun = @(X,Y,G) (separability_pcaLDA(X,Y,G,16,0), []); % wrap to return AUC if you refactor
% (or adapt to return a scalar AUC from your chosen routine)


% Plain grouped 5-fold:
cv = grouped_kfold(groups(trainMask|valMask), 5, 0);

% Or grouped + stratified by labels:
% cv = grouped_stratified_kfold(groups(trainMask|valMask), labels(trainMask|valMask), 5, 0);

% Example fold loop:
for f = 1:numel(cv)
    tr = cv{f}.train; va = cv{f}.val;
    Xtr = XTrain2(tr); Ytr = YTrain(tr);
    Xva = XTrain2(va); Yva = YTrain(va);
    % ...fit features on Xtr, train LDA/SVM, evaluate on Xva...
end

%%
% assume you already have: data (cell), labels (categorical), groups (cellstr)
% and a grouped split (trainMask, valMask), plus fs.

subjIDs = unique(groups,'stable');
models = struct;

for s = 1:numel(subjIDs)
    sid = subjIDs{s};
    idxTrain = trainMask & strcmp(groups, sid);
    idxVal   = valMask   & strcmp(groups, sid);

    if ~any(idxTrain) || ~any(idxVal), continue; end

    Xtr = data(idxTrain); Ytr = labels(idxTrain);
    Xva = data(idxVal);   Yva = labels(idxVal);

    % build template from this subjects TRAIN spikes only
    T = build_subject_template(Xtr(Ytr=='Spike'), fs, 'trim_ms', 250);

    % compute features
    Ftr = cell2mat(cellfun(@(x) feat_vs_template(single(x), T), Xtr, 'uni', false));
    Fva = cell2mat(cellfun(@(x) feat_vs_template(single(x), T), Xva, 'uni', false));

    % class weights to handle imbalance
    cats = categories(Ytr); cnt = countcats(Ytr);
    w = sum(cnt)./max(cnt,1); w = w/mean(w);
    classW = arrayfun(@(y) w(strcmp(cats, string(y))), Ytr);

    % train a small classifier (logistic; SVM also fine)
    mdl = fitclinear(Ftr, Ytr, ...
        'Learner','logistic', 'ClassNames',cats, ...
        'Weights', double(classW), 'Solver','lbfgs');

    % evaluate
    [Yp, scores] = predict(mdl, Fva);
    cm  = confusionmat(Yva, Yp, 'order', cats);
    acc = sum(diag(cm))/sum(cm,'all');
    rec = diag(cm)./max(1,sum(cm,2));
    bal = mean(rec);

    % store
    models.(sid) = struct('template',T, 'mdl',mdl, 'cm',cm, 'acc',acc, 'bal',bal);
    fprintf('%s  acc=%.2f  BA=%.2f\n', sid, acc, bal);
end

%%

grp = cellstr(string(groups));
subjIDs = unique(grp,'stable');
models = struct;

valFracSub = 0.2;   % per-subject validation fraction
seedSub    = 0;

for s = 1:numel(subjIDs)
    sid = subjIDs{s};
    idxSub = strcmp(grp, sid);
    if ~any(idxSub), continue; end

    Xs = data(idxSub);
    Ys = labels(idxSub);   % categorical: {'NoSpike','Spike'}

    % Subject-level split (ignores global masks)
    [trS, vaS] = subject_stratified_split(Ys, valFracSub, seedSub);

    if ~any(trS) || ~any(vaS)
        fprintf('[%s] Could not make both TRAIN and VAL; skipping.\n', sid);
        continue;
    end

    Xtr = Xs(trS);  Ytr = Ys(trS);
    Xva = Xs(vaS);  Yva = Ys(vaS);

    % Need at least 1 Spike in TRAIN to build template
    if ~any(Ytr=='Spike')
        fprintf('[%s] No TRAIN spikes; skipping.\n', sid);
        continue;
    end

    % Build template from TRAIN spikes
    T = build_subject_template(Xtr(Ytr=='Spike'), fs, 'trim_ms', 250);

    % Features (robust build; skips failed epochs)
    [Ftr, keepTr] = mkfeatmat(Xtr, T);
    [Fva, keepVa] = mkfeatmat(Xva, T);
    Ytr = Ytr(keepTr);  Yva = Yva(keepVa);

    if isempty(Ftr) || isempty(Fva)
        fprintf('[%s] Empty features after extraction; skipping.\n', sid);
        continue;
    end

    % Keep consistent class order
    cats = categories(Ytr);
    Ytr  = reordercats(Ytr, cats);
    Yva  = reordercats(Yva, cats);

    % Standardize using TRAIN only
    [Ftrz, muF, sigF] = sanitize_features(Ftr);
    Fvaz = sanitize_features(Fva, muF, sigF);

    % Class weights
    cnt = countcats(Ytr); w = sum(cnt)./max(cnt,1); w = w/mean(w);
    Wtr = arrayfun(@(y) w(strcmp(cats, string(y))), Ytr);

    % Train & eval
    mdl = fitclinear(Ftrz, Ytr, 'Learner','logistic', 'Solver','lbfgs', ...
                     'ClassNames', cats, 'Weights', double(Wtr));

    Yp = predict(mdl, Fvaz);


    [YvaU, YpU, cats] = unify_cats(Yva, Yp, categories(Ytr));
    cm  = confusionmat(YvaU, YpU, 'Order', cats);
    acc = sum(diag(cm))/sum(cm,'all');
    rec = diag(cm)./max(1,sum(cm,2)); bal = mean(rec);

    models.(sid) = struct('template',T,'mu',muF,'sig',sigF, ...
                          'mdl',mdl,'cm',cm,'acc',acc,'bal',bal, ...
                          'nTrain',sum(trS),'nVal',sum(vaS));
    fprintf('[%s] acc=%.2f  BA=%.2f  (TRAIN=%d VAL=%d)\n', sid, acc, bal, sum(trS), sum(vaS));
end





