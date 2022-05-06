%% Compute sensor level single trial power spectra
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = [f-1 f+1];
cfg.tapsmofrq    = 1;
cfg.keeptrials   = 'yes';
datapow           = ft_freqanalysis(cfg, cln_data);

%% identify the indices of trials with high and low alpha power
freqind = nearest(datapow.freq, f);
tmp     = datapow.powspctrm(:,:,freqind);
chanind = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));

%% compute the power spectrum for the median splitted data
cfg              = [];
cfg.trials       = indlow;
datapow_low      = ft_freqdescriptives(cfg, datapow);

cfg.trials       = indhigh;
datapow_high     = ft_freqdescriptives(cfg, datapow);

%% compute the difference between high and low
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'divide';
powratio      = ft_math(cfg, datapow_high, datapow_low);

%% plot the topography of the difference along with the spectra
cfg        = [];
cfg.layout = lay;
cfg.xlim   = [f-1 f+1];
figure; ft_topoplotER(cfg, powratio);

cfg         = [];
cfg.channel = {'MEG1342'};
ft_singleplotER(cfg, datapow_high, datapow_low);

%% Compute fourier spectra for frequency of interest according to the trial split
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 2;
cfg.foi        = f;

cfg.trials = indlow;
freq_low   = ft_freqanalysis(cfg, cln_data);

cfg.trials = indhigh;
freq_high  = ft_freqanalysis(cfg, cln_data);

%% compute sensor level Fourier spectra, to be used for cross-spectral density computation.
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = f;
freq           = ft_freqanalysis(cfg, cln_data);

%% compute the beamformer filters based on the entire data
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.grid              = individual_grid;
cfg.headmodel         = individual_headmodel;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '100%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.fixedori      = 'yes';
source = ft_sourceanalysis(cfg, freq);

% use the precomputed filters
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.grid              = individual_grid;
cfg.sourcemodel.filter       = source.avg.filter;
cfg.headmodel         = individual_headmodel;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
%         cfg.pcc.projectnoise  = 'yes';
cfg.pcc.fixedori      = 'yes';
source_low = ft_sourceanalysis(cfg, freq_low);
source_high = ft_sourceanalysis(cfg, freq_high);

source_low_avg  = ft_sourcedescriptives([], source_low);
source_high_avg = ft_sourcedescriptives([], source_high);

cfg           = [];
cfg.operation = 'log10(x1)-log10(x2)';
cfg.parameter = 'pow';
source_ratio  = ft_math(cfg, source_high_avg, source_low_avg);
source_ratio.pow(source_ratio.pow>0)=0;

%%
% create a fancy mask
% source_ratio.mask = (1+tanh(2.*(source_ratio.pow./max(source_ratio.pow(:))-0.5)))./2;
% 
% cfg = [];
% cfg.mask = 'pow';
% cfg.loc = 'max';
% cfg.template = mri_realigned;
% cfg.savefile = [];
% cfg.volnorm     = 2; % yes: 1
% cfg.method  = 'slice';
% vy_source_plot(cfg, source_ratio);

%% surface mapping
% s_in = source_ratio; mask = 'pow'; thre = 0.3;
% Run_surfacemap_template

%% Compute connectivity and network
% pause,
% close all,

cfg         = [];
cfg.method  ='coh';
cfg.complex = 'absimag';
source_conn_low = ft_connectivityanalysis(cfg, source_low);
source_conn_high = ft_connectivityanalysis(cfg, source_high);

cfg           = [];
cfg.operation = 'log10(x1)-log10(x2)';
cfg.parameter = 'cohspctrm';
conn_ratio  = ft_math(cfg, source_conn_high, source_conn_low);

cfg           = [];
cfg.method    = 'degrees';
cfg.parameter = 'cohspctrm';
cfg.threshold = .5;
net_ratio = ft_networkanalysis(cfg,conn_ratio);
net_ratio.mask = (1+tanh(2.*(net_ratio.degrees./max(net_ratio.degrees(:))-0.5)))./2;
net_ratio.dim  = source_conn_high.dim;

cfg = [];
cfg.mask = 'degrees';
cfg.loc = 'max';
cfg.template = mri_realigned;
cfg.savefile = [];
cfg.volnorm     = 2; % yes: 1
cfg.method  = 'slice';
cfg.nslices = [10,60];
vy_source_plot(cfg, net_ratio);

s_in = net_ratio; mask = 'degrees'; thre = 0.3;
Run_surfacemap_template