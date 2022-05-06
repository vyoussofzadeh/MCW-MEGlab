function [cln_data,report_ic] = vy_ica_cleaning_light6(cfg_main, data_in)

% satis = 0;
disp('ica cleaning ...');
if exist(cfg_main.savefile, 'file') && (cfg_main.overwrite) == 2
    load(cfg_main.savefile)
else
    
    savepath = fullfile(cfg_main.savepath,'ica');
    if exist(savepath, 'file') == 0, mkdir(savepath), end
    
    n = cfg_main.n; % ICs
    
    cfg = [];
    cfg.lay = cfg_main.lay;
    cfg.subj = cfg_main.subj;
    cfg.run = cfg_main.run;
    cfg.n = n;
    cfg.savefig = 0;
    cfg.allpath = cfg_main.allpath;
    [comp, freq] = vy_ica3(cfg, data_in);
    %     title(savepath)
    cfg = [];
    cfg.updatesens = 'no';
    
    %% flatness freq spectrum
    for jj=1:size(freq.powspctrm,1)
        num=geomean(freq.powspctrm(jj,:));
        den=mean(freq.powspctrm(jj,:));
        spf(jj)=num/den;
    end
    [~, fidx] = sort(spf,'ascend');
    freqbtrl = spf(fidx);
    thresh = 1.25.*min(spf);
    freqbtrl = find(spf < thresh);
    
    mfreq = mean(freq.powspctrm,2);mfreq = mfreq/(max(mfreq));
    [~, fidx] = sort(mfreq,'ascend');
    freqbtrl = spf(fidx);
    thresh = 1.25.*min(mfreq);
    freqbtrl = find(mfreq < thresh);
    %     fidx = fidx(1:3);
    
    %     mfreq = mean(freq.powspctrm,2);mfreq = mfreq/(max(mfreq));
    %     mfreq = mfreq/(max(mfreq));
    %     thresh = 1.1.*min(mfreq);
    %     freqbtrl = find(mfreq < thresh);
    %     [~, fidx] = sort(mfreq(freqbtrl),'ascend');
    %     freqbtrl = freqbtrl(fidx);
    %
    clear bch_freq_label_disp bch_freq_label
    for i=1:length(freqbtrl)
        bch_freq_label_disp{i,:} = comp.label{freqbtrl(i)};
        bch_freq_label{i,:} = ['-',comp.label{freqbtrl(i)}];
    end
    
    %%
    cfg = [];
    cfg.pflag    = 2; % yes:1, No:2
    cfg.saveflag = 2; % yes:1, No:2
    cfg.savepath = [];
    cfg.latency  = [comp.time{1}(1),comp.time{1}(end)];%[-200,900]./1000;
    cfg.rejectpercentage = .90;
    [~,report_ic] = vy_artifactreject4(cfg, comp);
    save('ica/ICAreport.mat','report_ic');
    disp(report_ic)
    
    %%
    report_ic_all = unique([report_ic.bchan; bch_freq_label_disp]);
    
    %% Rejecting bad trials, identified by the ICA
    trials = find(~ismember(1:length(comp.trial),report_ic.btrl));
    cfg = [];
    cfg.trials = trials;
    data_in = ft_selectdata(cfg, data_in);
    
    %%
    %     cfg = [];
    %     cfg.viewmode = 'component';
    %     cfg.layout = cfg_main.lay;
    %     ft_databrowser(cfg, r_data);
    %     colormap(brewermap(256, '*RdYlBu'));
    %     set(gcf, 'Position', [600   600   700   500]);
    
    %% Rejecting bad ICAs and prjecting back into data space
    %     disp('=============================')
    %     disp('suggested')
    %     disp(report_ic.bchan)
    %     cfg = [];
    %     if cfg_main.select == 1
    %         bic = input(['Select bad ICs for ' cfg_main.subj,':']);
    %         cfg.component = comp.label(bic);
    %
    %     else
    %         cfg.component =  report.bchan;
    %     end
    %     cfg.updatesens = 'no';
    %     cfg.trials = trials;
    %     cln_data = ft_rejectcomponent(cfg, comp, data_in);
    %     close all
    %     if cfg_main.saveflag == 1
    %         save(cfg_main.savefile, 'cln_data', '-v7.3');
    %     end
    
    %% Rejecting bad ICAs and prjecting back into data space
    warning off
    disp(['Select bad ICs for ' cfg_main.subj,':'])
    disp('=============================')
    disp('high Kurtosis')
    disp(report_ic.bchan)
    disp('Spectralflatness')
    disp(bch_freq_label_disp)
    disp('selected for rejection')
    disp(report_ic_all)
%     pause,
    cfg = [];
    if cfg_main.select == 1
        bic = input(['Select bad ICs for ' cfg_main.subj,':']);
        cfg.component = comp.label(bic);
    else
        
        inask = input ('look good, yes=1, no =2?');
        switch inask
            case 1
                cfg.component =  report_ic_all;
            case 2
                bic = input(['Select bad ICs for ' cfg_main.subj,':']);
                cfg.component = comp.label(bic);
        end
    end
    cfg.updatesens = 'no';
    cfg.trials = trials;
    cln_data = ft_rejectcomponent(cfg, comp, data_in);
    close all
    if cfg_main.saveflag == 1
        save(cfg_main.savefile, 'cln_data', '-v7.3');
        textfile_rej = 'ica/selected_badICs';
        badICs = cell2table(cfg.component);
        if ~isempty(badICs)
            badICs.Properties.VariableNames{'Var1'} = 'bICAs';
            %             writetable(badICs,textfile_rej,'Delimiter',' ');
        end
    end
    
    %% back to new ft
    restoredefaultpath
    addpath((cfg_main.allpath.ft_path));
    ft_defaults
    addpath(genpath(cfg_main.allpath.hcp_path));
    addpath(genpath(cfg_main.allpath.cd_org));
    addpath(genpath(cfg_main.allpath.exfig_path));
end

