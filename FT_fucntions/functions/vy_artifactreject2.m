function r_data = vy_artifactreject2(cfg_main, dat)

% disp('Identifying bad channels/sensors and bad trials ...');
% if exist(cfg_main.savepath, 'file') == 2
%     load(cfg_main.savepath)
% else

%% kurtosis
cfg = [];
cfg.trials = 'all';
cfg.metric = 'kurtosis';
cfg.channel = 'all';
cfg.latency = cfg_main.latency;
[level,info] = vy_compute_metric(cfg,dat);
%     metric.kurt = level;
info.pflag = cfg_main.pflag;
[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = vy_plot_chantrl(info,level);

thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel

%% zvalue
cfg = [];
cfg.trials = 'all';
cfg.metric = 'zvalue';
cfg.channel = 'all';
cfg.latency = cfg_main.latency;
[level,info] = vy_compute_metric(cfg,dat);
%     metric.kurt = level;
info.pflag = cfg_main.pflag;
[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = vy_plot_chantrl(info,level);

thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel

%% Var
cfg = [];
cfg.trials = 'all';
cfg.metric = 'var';
cfg.channel = 'all';
cfg.latency = cfg_main.latency;
[level,info] = vy_compute_metric(cfg,dat);
%     metric.kurt = level;
info.pflag = cfg_main.pflag;
[maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = vy_plot_chantrl(info,level);

thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel

%%
%     btrl_all = unique([btrl.kurtosis,btrl.var,btrl.zvalue]);
bch_all = unique([bch.kurtosis;bch.var;bch.zvalue]);

disp('Bad channels:')
for i=1:length(bch_all)
    bch_all_label_disp{i,:} = dat.label{bch_all(i)};
    bch_all_label{i,:} = ['-',dat.label{bch_all(i)}];
end
disp(bch_all_label_disp);

%% Removing bad trials
if cfg_main.rbadtrl == 1
    btrl_all = unique([btrl.kurtosis,btrl.var,btrl.zvalue]);
    disp('Bad trials:')
    disp(btrl_all);
    cfg = [];
    cfg.trials = find(~ismember(1:length(dat.trial),btrl_all));
    dat = ft_selectdata(cfg, dat);
    report.btrl = btrl_all;
end

%% Removing bad sensors
if cfg_main.rbadsen == 1
    if length(bch_all) < 10
        cfg = [];
        cfg.channel = ['all';bch_all_label];
        dat = ft_selectdata(cfg, dat);
    else
        warning('too many bad sensors w> 10, rejection skipped');
    end
    report.bchan = bch_all_label_disp;
end


r_data = dat;
if cfg_main.saveflag ==1
    %         save(cfg_main.savepath, 'r_data', 'report','-v7.3');
end
% end




%%
% cfg = [];
% cfg.trials = 'all';
% cfg.metric = 'kurtosis';
% cfg.channel = 'all';
% cfg.latency = [-400,900];
% [level,info] = vy_compute_metric(cfg,f_data2); metric.kurt = level;
% pflag = 1;
% [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = vy_plot_chantrl(info,level,pflag);
