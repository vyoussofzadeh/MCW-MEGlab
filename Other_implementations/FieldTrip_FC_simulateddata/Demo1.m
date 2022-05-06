
ft_path = '/data/MEG/Vahab/Github/fieldtrip';
addpath(ft_path);
ft_defaults

simcfg             = [];
simcfg.ntrials     = 500;
simcfg.triallength = 1;
simcfg.fsample     = 200;
simcfg.nsignal     = 3;
simcfg.method      = 'ar';
mu                 = 0.5;

% params(i,j,k): j -> i at t=k
simcfg.params(:,:,1) = [   0      0      0;
    1.0      0      0;
    0    1.0     mu];absnoise           = [ 1.0   0.2   0.3 ];


% paper defines stds, not cov:
simcfg.noisecov      = diag(absnoise.^2);
data2           = ft_connectivitysimulation(simcfg);

freq                   = [];
freq.freqcfg           = [];
freq.freqcfg.method    = 'mtmfft';
freq.freqcfg.output    = 'fourier';
freq.freqcfg.tapsmofrq = 2;
freqdata2           = ft_freqanalysis(freq.freqcfg, data2);

%%

cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.method      = 'ar';

cfg.params(:,:,1) = [ 0.8    0    0 ; 
                        0  0.9  0.5 ;
                      0.4    0  0.5];
                      
cfg.params(:,:,2) = [-0.5    0    0 ; 
                        0 -0.8    0 ; 
                        0    0 -0.2];
                        
cfg.noisecov      = [ 0.3    0    0 ;
                        0    1    0 ;
                        0    0  0.2];

data              = ft_connectivitysimulation(cfg);


%%
% figure
% plot(data.time{1}, data.trial{1}) 
% legend(data.label)
% xlabel('time (s)')

%%
cfg = [];
cfg.viewmode = 'vertical';  % you can also specify 'butterfly' 
ft_databrowser(cfg, data);

%%
cfg         = [];
cfg.order   = 5;
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, data);

%%
cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);

%%
cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.tapsmofrq = 2;
freq          = ft_freqanalysis(cfg, data);

%%
cfg           = [];
cfg.method    = 'coh';
coh           = ft_connectivityanalysis(cfg, freq);
cohm          = ft_connectivityanalysis(cfg, mfreq);

cfg           = [];
cfg.parameter = 'cohspctrm';
cfg.zlim      = [0 1];
ft_connectivityplot(cfg, coh, cohm);

%%
cfg           = [];
cfg.method    = 'plv';
plv           = ft_connectivityanalysis(cfg, freq);
plvm          = ft_connectivityanalysis(cfg, mfreq);

figure
cfg           = [];
cfg.parameter = 'plvspctrm';
cfg.zlim      = [0 1];
ft_connectivityplot(cfg, plv, plvm);

%%
cfg           = [];
cfg.method    = 'wpli_debiased';
wpli           = ft_connectivityanalysis(cfg, freq);
wplim          = ft_connectivityanalysis(cfg, mfreq);
figure
cfg           = [];
cfg.parameter = 'wpli_debiasedspctrm';
cfg.zlim      = [0 1];
ft_connectivityplot(cfg, wpli, wplim);

%%
cfg           = [];
cfg.method    = 'pdc';
pdc           = ft_connectivityanalysis(cfg, freq);
pdcm          = ft_connectivityanalysis(cfg, mfreq);

figure
cfg           = [];
cfg.parameter = 'pdcspctrm';
cfg.zlim      = [0 1];
ft_connectivityplot(cfg, pdc, pdcm);

%%
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, mfreq);

figure,
cfg           = [];
cfg.parameter = 'grangerspctrm';
cfg.zlim      = [0 1];
ft_connectivityplot(cfg, granger);