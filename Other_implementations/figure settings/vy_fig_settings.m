
%% Matlab figure settings
figure,
bar((DL_runs.math + DL_runs.str)/2), L = length(DL_runs.math);
set(gca,'Xtick', 1:L,'XtickLabel',sub_runs_name);
set(gca,'FontSize',10,'XTickLabelRotation',90);
set(gcf, 'Position', [1000   100   1500   300]);
set(gca,'color','none');
title([tag, ', DL, mean (Str, Math) '])

%% Legend
legend('hide')

%% Axis labelling
xlabel(rois(idx(1)))
ylabel([tag,'-',ttag])

%% Corrolation
[R,P] = corrcoef(X(:,idx(1)),y1');

%% Regression
tbl = table(X(:,idx(1)),y1');
mdl = fitlm(tbl,'linear');
plot(mdl)

%%
