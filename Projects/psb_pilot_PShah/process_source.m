%% The psb_pilot

% MEG phase amplitude coupling
% Writtern by MCW group, Shah-Basak, Priyanka <prishah@mcw.edu>
% input data from MEG (pre-) processing pipeline by Youssofzadeh, Vahab
% Lastest update: 05/25/2022

clear; clc, close('all'); warning off

%% Paths
restoredefaultpath
maindir = '/group/prishah/LanguageMEG/psb_pilot';
ssid = 'ss_pilot_2';
subjid = 'pilot2';
datcol = '220512';

script_path = maindir;%'/data/MEG/Research/psb_pilot';
addpath(genpath(script_path));
%- Input dir
indir = fullfile(maindir, ssid, datcol, 'tsss');%'/data/MEG/Research/psb_pilot/ss_pilot_2/220512/tsss';
inpath = fullfile(indir, 'STM_Block1_raw.fif');

%- Output dir
outdir = fullfile(maindir, 'FT');
%- MRI dir
mridir = fullfile(maindir, ssid, datcol, 'mri');
%
ft_path = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_20190419');
addpath(ft_path);
ft_defaults

allpath = [];
allpath.ft_path18 = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_041718'); % needed for IC plotting
allpath.ft_path = ft_path;
allpath.script_path = script_path;

%% Load data
taskd = 'tone'; %options: 'pstm' or 'tone'
load(fullfile(outdir, ['all_dataclean_', taskd, '_', subjid, '.mat'])) %loads alldata
ftdatatone = alldata;

taskd = 'pstm'; %options: 'pstm' or 'tone'
load(fullfile(outdir, ['all_dataclean_', taskd, '_', subjid, '.mat'])) %loads alldata
ftdatastm = alldata;
%% Cross-spectral density
% cfg = [];
% cfg.keepsampleinfo = 'no';
% dataBoth = ft_appenddata(cfg,ftdatastm, ftdatatone);
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'powandcsd';
cfg.tapsmofrq = 4;
cfg.foilim = [8 8];
freqstm = ft_freqanalysis(cfg, ftdatastm);
freqtone = ft_freqanalysis(cfg, ftdatatone);
%freqBoth = ft_freqanalysis(cfg,dataBoth);
%% Read Anatomy and segment the brain
mri = ft_read_mri(fullfile(mridir,'mprage.nii.gz'));
mri.coordsys = 'ras';
%mri_coordsys = ft_determine_coordsys(mri);
ft_sourceplot([],mri)
%% Volume realign based on fiducials
cfg = [];
cfg.method   = 'interactive';
cfg.coordsys = 'neuromag';
mri_realigned = ft_volumerealign(cfg, mri);
%% Read headshape
inpath = fullfile(indir, 'STM_Block1_raw.fif');
headshape = ft_read_headshape(inpath)
figure; ft_plot_headshape(headshape)
%% Coregister with digitized head points
grad    = ft_read_sens(inpath,'senstype','meg'); % Load MEG sensors
%elec    = ft_read_sens(inpath,'senstype','eeg'); % Load EEF electrodes

% Plot head points
figure;
ft_plot_headshape(headshape);     % Plot headshape again
ft_plot_sens(grad);
%ft_plot_sens(elec, 'style', '*b');
%ft_sourceplot([],mri)
 
%% Volume realign method 2

cfg = [];
cfg.method              = 'headshape';
cfg.headshape.headshape = headshape;
cfg.headshape.icp       = 'yes';
cfg.coordsys            = 'neuromag';

mri_realigned_2 = ft_volumerealign(cfg, mri_realigned); %_realigned

%% checking co-registration
cfg = [];
cfg.method              = 'headshape';
cfg.headshape.headshape = headshape;
cfg.coordsys            = 'neuromag';
cfg.headshape.icp       = 'no';        %Do not fit point again
mri_realigned_3 = ft_volumerealign(cfg, mri_realigned_2);
%% Reslice MRI data 
cfg = [];
cfg.resolution = 1;
mri_resliced = ft_volumereslice(cfg, mri_realigned); %mri_realigned_3
mri_resliced = ft_convert_units(mri_resliced, 'cm');
%% Segmenting the brain1
disp('segmenting the brain');
cfg = [];
cfg.output = {'brain','skull','scalp'};
cfg.write = 'no';
%cfg.name = 'brain';
segmentedmri = ft_volumesegment(cfg,mri_resliced);

% Check the segmentation
cfg = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg, segmentedmri);
cfg.funparameter = 'skull';
ft_sourceplot(cfg, segmentedmri);
cfg.funparameter = 'scalp';
ft_sourceplot(cfg, segmentedmri);
%% Construct meshes
cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = 'brain';
cfg.numvertices = 3000;

mesh_brain = ft_prepare_mesh(cfg, segmentedmri);

%%
disp('generating the head model')
cfg = [];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, mesh_brain);


%% Visualize
figure 
ft_plot_sens(freqstm.grad,'unit','cm');
hold on
ft_plot_headshape(headshape, 'unit', 'cm')
ft_plot_headmodel(headmodel, 'unit', 'cm')
ft_plot_axes([], 'unit', 'cm');
%% Source model and lead fields
cfg= [];
cfg.grad = freqstm.grad;
cfg.headmodel = headmodel;
cfg.sourcemodel.unit = 'cm';
cfg.normalize = 'yes';
sourcemodel = ft_prepare_leadfield(cfg);
%% Compute inverse filters
cfg = []; 
cfg.method = 'dics';
cfg.frequency = 8;
cfg.sourcemodel = sourcemodel;
cfg.headmodel = headmodel;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda = '5%';
cfg.dics.keepfilter = 'yes';
cfg.dics.realfiler = 'yes';
sourcestm = ft_sourceanalysis(cfg, freqstm);
sourcetone = ft_sourceanalysis(cfg, freqtone);
sourceDiff = sourcestm;
sourceDiff.avg.pow = (sourcestm.avg.pow ./ sourcestm.avg.noise) - (sourcetone.avg.pow ./ sourcetone.avg.noise);
sourceNAI = sourcestm;
sourceNAI.avg.pow = sourcestm.avg.pow ./ sourcestm.avg.noise;
%% Interpolate MRI 
cfg = [];
cfg.downsample = 2;
cfg.parameter= 'pow';
sourceInt = ft_sourceinterpolate(cfg,sourceNAI,mri_resliced);
%% Visualize
maxval = max(sourceInt.pow)
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim = [0 maxval];
cfg.opacitylim = [0 maxval];
cfg.opacitymap = 'rampup';
ft_sourceplot(cfg, sourceInt)
%% convert to a standard template
cfg = [];
cfg.nonlinear = 'yes';
sourceIntNorm = ft_volumenormalise(cfg,sourceInt);
%%
maxval = max(sourceInt.pow)
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim = [0 maxval];
cfg.opacitylim = [0 maxval];
cfg.opacitymap = 'rampup';
ft_sourceplot(cfg, sourceIntNorm)
%%
cfg = [];
cfg.method = 'surface';
cfg.funparameter = 'pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim = [0 maxval];
cfg.funcolormap = 'hsv';
cfg.opacitylim = [0 maxval];
cfg.opacitymap = 'rampup';
cfg.projmethod = 'nearest';
cfg.surffile = 'surface_white_both.mat';
cfg.surfdownsample = 5;
cfg.interactive = 'yes';
ft_sourceplot(cfg,sourceIntNorm);
view([90 0])