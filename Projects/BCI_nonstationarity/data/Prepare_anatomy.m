%% BCI dataset, Ulster University & Medical College of Wisconsin

% Script: Anatomy prepare
% Project: BCI_nonstationarity
% Writtern by: Vahab Youssof Zadeh
% Update: 11/01/2022


load('standard_mri')
load temp_grid
load('standard_singleshell')

cfg                 = [];
cfg.grid.warpmni    = 'yes';
cfg.spmversion     = 'SPM12';
cfg.grid.nonlinear  = 'yes';
cfg.grid.template   = template_grid;
cfg.mri             = mri;
cfg.grid.unit       = 'mm';
individual_grid     = ft_prepare_sourcemodel(cfg);


vol = ft_convert_units(vol,'mm');
individual_headmodel = vol;

figure;
ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
hold on;
ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
view ([0 90])

save('anat.mat', 'individual_grid', 'individual_headmodel')