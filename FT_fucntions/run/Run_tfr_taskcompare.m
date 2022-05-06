%%
tsk = 'DFN';
load(fullfile(tsk,[tsk,'_tfr.mat']))
tfr_dfn = tfr_all;
sub_dfn = Sub_all;

tsk = 'PN';
load(fullfile(tsk,[tsk,'_tfr.mat']))
tfr_pn = tfr_all;
sub_pn = Sub_all;
clear tfr_all Sub_all names

% The finer time and frequency axes:
tim_interp = linspace(-0.5, 2, 512);
freq_interp = linspace(1, 40, 512);
[tim_grid_interp, freq_grid_interp] = meshgrid(tim_interp, freq_interp);

%% Group, DFN
pow_interp_dfn_all = []; time_of_interest = []; freq_of_interest = [];
nsub = length(tfr_dfn);
for i=1:nsub

    [tim_grid_orig, freq_grid_orig]     = meshgrid(tfr.time, linspace(1,40, size(tfr_dfn{i},1)));
    pow_interp = interp2(tim_grid_orig, freq_grid_orig, tfr_dfn{i}, tim_grid_interp, freq_grid_interp, 'spline');
    pow_interp_dfn_all(i,:,:) = (pow_interp);
    
    pow_interp1  = pow_interp(50:end,50:end);
    tim_interp1  = tim_interp(50:end);
    freq_interp1 = freq_interp(50:end);    
    [~,idx] = min(pow_interp1(:));
    [row,col] = ind2sub(size(pow_interp1),idx);   
    time_of_interest(i) = tim_interp1(col);
    freq_of_interest(i) = freq_interp1(row);
    
end
pow_interp = squeeze(mean(pow_interp_dfn_all,1));

toi_dfn = [num2str(mean(time_of_interest)),'+-', num2str(std(time_of_interest))];
disp(toi_dfn)

tfr_sel = pow_interp;
% [h,p,ci,stats] = ttest(pow_interp)
Run_plot_tfr

%%
load('/data/MEG/Clinical/group_dics/taskcomp/baddata.mat');

%%
[C,ia,ib] = intersect(C1,sub_dfn, 'stable');
pow_interp = squeeze(mean(pow_interp_all(ib,:,:),1));
Run_plot_tfr

%% Group, PN
pow_interp_pn_all = []; time_of_interest = []; freq_of_interest = [];
for i=1:length(tfr_pn)
    [tim_grid_orig, freq_grid_orig]     = meshgrid(tfr.time, linspace(1,40, size(tfr_pn{i},1)));
    pow_interp = interp2(tim_grid_orig, freq_grid_orig, tfr_pn{i}, tim_grid_interp, freq_grid_interp, 'spline');
    pow_interp_pn_all(i,:,:) = pow_interp;
%     pow_interp_all(i,:,:) = pow_interp./max(pow_interp(:));

    pow_interp1  = pow_interp(50:end,50:end);
    tim_interp1  = tim_interp(50:end);
    freq_interp1 = freq_interp(50:end);    
    [~,idx] = min(pow_interp1(:));
    [row,col] = ind2sub(size(pow_interp1),idx);   
    time_of_interest(i) = tim_interp1(col);
    freq_of_interest(i) = freq_interp1(row);
end
pow_interp = squeeze(mean(pow_interp_pn_all,1));

toi_pn = [num2str(mean(time_of_interest)),'+-', num2str(std(time_of_interest))]
disp(toi_pn)

Run_plot_tfr

%%
[C,ia,ib] = intersect(C1,sub_pn, 'stable');
pow_interp = squeeze(mean(pow_interp_pn_all(ib,:,:),1));
Run_plot_tfr

%% TFR subject
sel = 1;
tfr_sel = tfr_dfn{sel};
disp(sub_dfn{1})

[tim_grid_orig, freq_grid_orig]     = meshgrid(tfr.time, linspace(1,40, size(tfr_sel,1)));
pow_interp = interp2(tim_grid_orig, freq_grid_orig, tfr_sel, tim_grid_interp, freq_grid_interp, 'spline');

Run_plot_tfr

%% Stats,

% baseline=mean(TFrAll.powspctrm(:,:,:,1:6),4);
% for timei=2:31;
%     TFrAll.powspctrm(:,:,:,timei)=TFrAll.powspctrm(:,:,:,timei)-baseline;
% end
% % no compute the statistic
% cfg = [];
% cfg.method = 'stats';
% cfg.design(1,:) = [ones(1,nsub)];
% cfg.latency     = [0 0.35];
% cfg.frequency   = [1 20];
% cfg.statistic = 'ttest'; % compares the mean to zero
% cfg.feedback='no';
% frstat=ft_freqstatistics(cfg,TFrAll);
% % now plot 1-probability (1 = sig, less than 0.95 not sig)
% cfg=[];
% cfg.layout = '4D248.lay';
% frstat.powspctrm = 1-frstat.prob;
% cfg.zlim=[0.999 1]
% cfg.interactive='yes';
% fig3=figure;
% set(fig3,'Position',[0,0,800,800]);
% ft_multiplotTFR(cfg, frstat);



