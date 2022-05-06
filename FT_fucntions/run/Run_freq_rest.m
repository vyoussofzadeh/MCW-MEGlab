%% freq analysis (fft)

savepath = fullfile(outd.sub,'Freq');
if exist(savepath, 'file') == 0, mkdir(savepath), end

%%
clc,
%- checking inter-trial time-intervals
TI = datain.time{2}(1) - datain.time{1}(1);
datain.time{1}(end) - datain.time{1}(1);

if exist(fullfile(savepath,'tfr.mat'),'file')~= 2
    
    disp('Enter the highest freq in the data in Hz, e.g. 40:')
%     fmax = input(':');
    fmax = 40;
    % do tfr-decomposition
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    % cfg.taper      = 'dpss';
    cfg.foi        = 1:2:fmax;
    cfg.keeptrials = 'yes';
    cfg.t_ftimwin  = 3./cfg.foi;
    % cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.tapsmofrq  = 0.8 *cfg.foi;
    cfg.toi        = datain.time{1}(1):0.05:datain.time{1}(end);
    tfr            = ft_freqanalysis(cfg, datain);
    save(fullfile(savepath,'tfr.mat'),'tfr','fmax');
    if TI ~= 0, tfr.time = linspace(0, datain.time{1}(end) - datain.time{1}(1), length(cfg.toi)); end
    
else
    load(fullfile(savepath,'tfr.mat'));
    if ~exist('fmax','var'), fmax = 30; end
end

tfr.powspctrm(isnan(tfr.powspctrm))=0;
tfr.time = linspace(0, tfr.time(end) - tfr.time(1), length(tfr.time));


cfg = [];
cfg.savepath = 1;
cfg.savefile = fullfile(savepath,'tfr');
cfg.fmax = fmax;
cfg.toi = [tfr.time(1), tfr.time(end)];
cfg.bslcorr = 2;
cfg.effect = 1; % postive = 1, negative = 2; both = 3;
[time_of_interest,freq_of_interest] = vy_tfr_plot_rest(cfg, tfr);
disp(['time_of_interest:',num2str(time_of_interest),'sec']);
disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
L = 0.3;

