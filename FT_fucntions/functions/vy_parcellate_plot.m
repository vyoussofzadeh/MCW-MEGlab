function sourceint_pow = vy_parcellate_plot(data_intpar,coor, mask)

%- surface visualisation
% cfg               = [];
% cfg.method        = 'surface';
% cfg.funparameter  = 'anatomy';
% cfg.maskparameter = cfg.funparameter;
% cfg.funcolorlim   = [-0.3 0.3];
% cfg.opacitylim    = [-8 8];
% cfg.opacitymap    = 'rampup';
% cfg.funcolormap   = 'jet';
% cfg.colorbar      = 'yes';
% cfg.projthresh     = 0.6;
% ft_sourceplot(cfg, data_intpar);
% view ([-100 0 0]), light ('Position',[-100 0 0])

%%
tmp = abs(data_intpar.(mask));
tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
data_intpar.(mask) = tmp;

%%
sMRI = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
cfg1 = [];
cfg1.sourceunits  = 'mm';
cfg1.parameter    = mask;
cfg1.downsample   = 1;
% cfg.interpmethod = 'sphere_avg';
% cfg1.coordsys     = 'mni';
sourceint_pow = ft_sourceinterpolate(cfg1, data_intpar, ft_read_mri(sMRI, 'format', 'nifti_spm'));

% sourceint_pow.anatomy(isnan(sourceint_pow.anatomy(:))) = 0;
sourceint_pow.(mask)(isnan(sourceint_pow.(mask)(:))) = 0;

%%

clear savepath
savepath{1} = [mask,'_left'];
savepath{2} = [mask,'_right'];
cfg = [];
cfg.subj = 'par';
cfg.mask = mask;
cfg.thre = 0.70;
cfg.savepath = savepath;
cfg.saveflag = 2;
vy_mapvisualisation(cfg, sourceint_pow);
% vy_mapvisualisation(sourceint_pow,'anatomy',0.7, savepath);

vy_ROI_report(data_intpar,.75, coor, mask);
savepath = ['ROI_',mask];
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300)

