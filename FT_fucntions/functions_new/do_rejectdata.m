function [r_data,report] = do_rejectdata(cfg_main, dat)

disp('Bad channel/trial ...');
if exist(cfg_main.savepath, 'file') == 2
    load(cfg_main.savepath)
else
    switch cfg_main.method
        case 'auto'
            %% kurtosis
            cfg = [];
            cfg.trials = 'all';
            cfg.metric = 'kurtosis';
            cfg.channel = 'all';
            cfg.latency = cfg_main.latency;
            [level,info] = do_compute_metric(cfg,dat);
            %     metric.kurt = level;
            info.pflag = cfg_main.pflag;
            [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = do_plot_chantrl(info,level);
            
            thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
            thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel
            
            %% zvalue
            cfg = [];
            cfg.trials = 'all';
            cfg.metric = 'zvalue';
            cfg.channel = 'all';
            cfg.latency = cfg_main.latency;
            [level,info] = do_compute_metric(cfg,dat);
            %     metric.kurt = level;
            info.pflag = cfg_main.pflag;
            [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = do_plot_chantrl(info,level);
            
            thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
            thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel
            
            %% Var
            cfg = [];
            cfg.trials = 'all';
            cfg.metric = 'var';
            cfg.channel = 'all';
            cfg.latency = cfg_main.latency;
            [level,info] = do_compute_metric(cfg,dat);
            %     metric.kurt = level;
            info.pflag = cfg_main.pflag;
            [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = do_plot_chantrl(info,level);
            
            thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxpertrl_all); btrl.(cfg.metric) = find(maxpertrl > thresh.(cfg.metric)); % Trials
            thresh.(cfg.metric) = cfg_main.rejectpercentage.*max(maxperchan_all); bch.(cfg.metric) = find(maxperchan > thresh.(cfg.metric)); % Channel
            
            %%
            btrl_all = unique([btrl.kurtosis,btrl.var,btrl.zvalue]);
            %     btrl_all = unique([btrl.kurtosis,btrl.var,btrl.zvalue]);
            bch_all = unique([bch.kurtosis;bch.var;bch.zvalue]);
            disp('Bad trials:')
            disp(btrl_all);
            disp('Bad channels:')
            for i=1:length(bch_all)
                bch_all_label_disp{i,:} = dat.label{bch_all(i)};
                bch_all_label{i,:} = ['-',dat.label{bch_all(i)}];
            end
            disp(bch_all_label_disp);
            
            %% Removing bad channels/trials
            cfg = [];
            cfg.trials = find(~ismember(1:length(dat.trial),btrl_all));
            dat = ft_selectdata(cfg, dat);
            
            report.btrl = btrl_all; report.bchan = bch_all_label_disp;
            r_data = dat;
            
            if ~isempty(report)
                disp('REJECTED TRIALS:')
                disp(report.btrl')
            end
        case 'manual'
            cfg = [];
            cfg.metric = 'kurtosis';  % use by default kurtosis method
            cfg.latency = [dat.time{1}(1),dat.time{1}(end)];
            r_data   = ft_rejectvisual(cfg, dat);
            
            % Bad trial info
            data = ft_checkdata(dat, 'datatype', {'raw+comp', 'raw'}, 'feedback', 'yes', 'hassampleinfo', 'yes');
            btrlsample = r_data.cfg.artfctdef.summary.artifact;
            if ~isempty(btrlsample)
                for l=1:size(btrlsample,1)
                    btrl(l,:) = find(ismember(data.sampleinfo(:,1),btrlsample(l))==1);
                end
                report.btrl = btrl;
                disp('REJECTED TRIALS:')
                disp(btrl)
            end
    end
    if cfg_main.saveflag ==1
        %         save(cfg_main.savepath, 'r_data', 'report','-v7.3');
    end
end


%%
% cfg = [];
% cfg.trials = 'all';
% cfg.metric = 'kurtosis';
% cfg.channel = 'all';
% cfg.latency = [-400,900];
% [level,info] = do_compute_metric(cfg,f_data2); metric.kurt = level;
% pflag = 1;
% [maxperchan, maxpertrl, maxperchan_all, maxpertrl_all] = do_plot_chantrl(info,level,pflag);
