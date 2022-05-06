%%
f = input('FOI (Hz)? ');
tapsmofrq = input('tapsmofrq, e.g. 4 Hz? ');

cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [f f];
cfg.plotflag  = 2;
cfg.taper    = 'dpss'; cfg.tapsmofrq  = tapsmofrq;

if f < 4, cfg.tapsmofrq  = 1; cfg.taper    = 'hanning'; end

f_data = [];
[f_data.pst,~,~,~] = vy_fft(cfg, cln_data); f_data.pst.elec = cln_data.grad;
cfg.foilim = [15 15];
cfg.tapsmofrq  = 12;
[f_data.bsl,~,~,~] = vy_fft(cfg, cln_data); f_data.bsl.elec = cln_data.grad;

%% DICS analysis (relative power)
cfg = [];
cfg.method = 'dics';
cfg.dics.lambda = '100%';
cfg.grid = individual_grid;
cfg.frequency    = f_data.pst.freq;
cfg.headmodel = individual_headmodel;
cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);

tmp = s_data.bsl.avg.pow; tmp(isnan(tmp))=0; tmp = tmp./max(tmp); s_data.bsl.avg.pow = tmp;
tmp = s_data.pst.avg.pow; tmp(isnan(tmp))=0; tmp = tmp./max(tmp); s_data.pst.avg.pow = tmp;

cfg = [];
cfg.parameter = 'avg.pow';
cfg.operation = '(x1-x2)';
source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
source_diff_dics.pow = (s_data.pst.avg.pow - s_data.bsl.avg.pow)./(s_data.pst.avg.pow + s_data.bsl.avg.pow);
source_diff_dics.pow(source_diff_dics.pow<0)=0;

cfg = [];
cfg.mask = 'pow';
cfg.loc = 'max';
cfg.template = mri_realigned;
cfg.savefile = [];
cfg.volnorm     = 2; % yes: 1
        cfg.method  = 'ortho';
% cfg.method  = 'slice';
cfg.nslices = [10,60];
vy_source_plot(cfg, source_diff_dics);


%% surface mapping
s_in = source_diff_dics; mask = 'pow'; thre = 0;
Run_surfacemap_template