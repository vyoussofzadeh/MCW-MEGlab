function stat_group = vy_betweensub_stat(pow,DFN2)

num_sub = size(pow,1);
% subject loop
for k = 1:num_sub
    stat_all{k}.pow = pow(k,:)';
    stat_all{k}.pos = DFN2.source_diff_dics.pos;
    stat_all{k}.inside = DFN2.source_diff_dics.inside;
    stat_all{k}.dim = DFN2.source_diff_dics.dim;
end

%
% Create the null structure
% data_N = source_diff_dics;
data_N = [];
for k = 1:num_sub
    data_N{k}.pow(:) = zeros(size(pow,2),1);
    data_N{k}.pos = DFN2.source_diff_dics.pos;
    data_N{k}.inside = DFN2.source_diff_dics.inside;
    data_N{k}.dim = DFN2.source_diff_dics.dim;
    %data_N{k}.stat(:) = nanmean(stat_all{k}.stat(:));
end

% specify design matrix
design                           = zeros(2, 2*num_sub);
design(1, 1:num_sub)             = 1:num_sub;
design(1, num_sub + 1:num_sub*2) = 1:num_sub;
design(2, 1:num_sub)             = 1;
design(2, num_sub + 1:num_sub*2) = 2;
% design                           = 1:num_sub;

% second-level t-test
cfg                  = [];
cfg.method           = 'montecarlo';
cfg.parameter        = 'pow';
cfg.statistic        = 'depsamplesT';
% cfg.statistic        = 'indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.alpha            = 0.05; % adjust alpha-level for two-sided test
cfg.correcttail      = 'prob';
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the
% cfg.tail        = 0;
cfg.numrandomization = 5000;
cfg.design           = design;
cfg.uvar             = 1; % the 1st row in cfg.design contains the subject number (or trial number)
cfg.ivar             = 2; % the 2nd row in cfg.design contains the independent variable
stat_group           = ft_sourcestatistics(cfg, stat_all{:}, data_N{:});


end