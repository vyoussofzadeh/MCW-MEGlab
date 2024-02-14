function [freq, ff, psd,tapsmofrq] = do_fft(cfg_main, data)

%
cfg              = [];
cfg.method       = 'mtmfft';
cfg.output       = 'fourier';
cfg.keeptrials   = 'yes';
cfg.foilim       = cfg_main.foilim;
cfg.tapsmofrq    = cfg_main.tapsmofrq;
cfg.taper        = cfg_main.taper; %'hanning';
cfg.pad          = cfg_main.pad; %4;
freq             = ft_freqanalysis(cfg, data);
psd = squeeze(mean(mean(abs(freq.fourierspctrm),2),1));
ff = linspace(1, cfg.foilim(2), length(psd));

if cfg_main.plotflag ==1
    figure,plot(ff,psd)
    xlabel('Hz'); ylabel('psd')
end

tapsmofrq = cfg.tapsmofrq;

if cfg_main.saveflag ==1
    hcp_write_figure([cfg_main.savefile,'.png'], gcf, 'resolution', 300); 
end

