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
[f_data.pst,~,~,~] = do_fft(cfg, ep_data.pst); f_data.pst.elec = data_clean.grad;
[f_data.bsl,~,~,~] = do_fft(cfg, ep_data.bsl); f_data.bsl.elec = data_clean.grad;

%- DICS analysis (relative power)
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
cfg.operation = '(x1-x2)/(x1+x2)';
source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
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
do_source_plot(cfg, source_diff_dics);


source_tmp = source_diff_dics;
source_tmp.pos     = template_grid.pos;
source_tmp.dim     = template_grid.dim;
source_tmp.inside  = template_grid.inside;

cfg = [];
cfg.parameter = 'pow';
% cfg.interpmethod = 'sphere_avg';
cfg.interpmethod = 'smudge';
cfg.coordsys     = 'mni';
source_int = ft_sourceinterpolate(cfg, source_tmp, template_mri);

cfg = [];
cfg.subj = [num2str(f), 'Hz'];
cfg.mask = 'pow';
cfg.thre = 0.3;
cfg.savepath = savepath;
cfg.colorbar = 2;
cfg.saveflag = 2;
cfg.colormap = brewermap(256, '*RdYlBu');
cfg.surfinflated   = 'surface_inflated_both.mat';
cfg.views = [90,0; -90, 0];
do_mapvis(cfg, source_int);