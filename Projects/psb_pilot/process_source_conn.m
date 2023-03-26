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
%outdir = fullfile(maindir, 'FT');
ftoutdir = fullfile(maindir, ssid, datcol, 'meg');
%- MRI dir
mridir = fullfile(maindir, ssid, datcol, 'mri');
ftanatdir = fullfile(maindir, ssid, datcol, 'meg' , 'anat');
%
ft_path = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip-20210517'); %fieldtrip_20190419
addpath(ft_path);
ft_defaults

%addpath('/opt/matlab_toolboxes/ft_packages/fieldtrip-20210517/')
%ft_defaults
        
allpath = [];
allpath.ft_path18 = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_041718'); % needed for IC plotting
allpath.ft_path = ft_path;
allpath.script_path = script_path;

%% Load data
taskd = 'tone'; %options: 'pstm' or 'tone'
load(fullfile(ftoutdir, ['all_dataclean_', taskd, '_', subjid, '.mat'])) %loads alldata
ftdatatone = alldata;

taskd = 'pstm'; %options: 'pstm' or 'tone'
load(fullfile(ftoutdir, ['all_dataclean_', taskd, '_', subjid, '.mat'])) %loads alldata
ftdatastm = alldata;

%% template grid
load(fullfile(ft_path, 'template/sourcemodel/standard_sourcemodel3d10mm'));
template_grid = sourcemodel;
clear sourcemodel;

%% Load Anatomy, resliced and neuromag anatomy, segmented brain, headshape
load(fullfile(ftanatdir, 'headshape.mat'))
load(fullfile(ftanatdir, 'mri_resliced.mat'))
load(fullfile(ftanatdir, 'segmentedmri.mat'))

%% Coregister with digitized head points
grad    = ft_read_sens(inpath,'senstype','meg'); % Load MEG sensors
%elec    = ft_read_sens(inpath,'senstype','eeg'); % Load EEF electrodes

% Plot head points
figure;
ft_plot_headshape(headshape);     % Plot headshape again
ft_plot_sens(grad);
%ft_plot_sens(elec, 'style', '*b');
%ft_sourceplot([],mri)
 
% Check the segmentation
cfg = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg, segmentedmri);

%%
disp('generating the head model')
cfg = [];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, segmentedmri);

%% Source model and lead fields
cfg= [];
cfg.warpmni = 'yes';
cfg.template = template_grid;
cfg.nonlinear = 'yes';
cfg.mri = mri_resliced;
cfg.unit = 'cm';
sourcemodel = ft_prepare_sourcemodel(cfg);

%% Visualize
figure 
ft_plot_sens(grad,'unit','cm');
hold on
ft_plot_headshape(headshape, 'unit', 'cm')
ft_plot_headmodel(headmodel, 'unit', 'cm','facecolor','cortex','edgecolor','none')
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
ft_plot_axes([], 'unit', 'cm');
%% Cross-spectral density
fq = 8;
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'powandcsd';
cfg.tapsmofrq = 5;
cfg.foilim = [fq fq];
cfg.keeptrials = 'yes';
freqstm = ft_freqanalysis(cfg, ftdatastm);
freqtone = ft_freqanalysis(cfg, ftdatatone);

%% DICS
cfg = []; 
cfg.method = 'dics';
cfg.frequency = fq;
cfg.sourcemodel = sourcemodel;
cfg.headmodel = headmodel;
cfg.dics.projectnoise = 'yes';
%cfg.dics.lambda = '10%';
cfg.dics.keepfilter = 'yes';
cfg.dics.realfiler = 'yes';
%cfg.rawtrial = 'yes';
%cfg.keeptrials = 'yes';
sourcestm = ft_sourceanalysis(cfg, freqstm);
%save(fullfile(ftoutdir,'sample_source_singletrial.mat'),'sourcestm','-v7.3');
sourcetone = ft_sourceanalysis(cfg, freqtone);
%%
cfg = [];
sourcestm =ft_sourcedescriptives(cfg, sourcestm);
sourcetone =ft_sourcedescriptives(cfg, sourcetone);
cfg = [];
cfg.parameter = 'pow';
cfg.operation = 'log(x1/x2)';
sourceDiff = ft_math(cfg, sourcestm, sourcetone);
%% Interpolate MRI 
cfg = [];
cfg.parameter= 'nai';
sourcestmInt = ft_sourceinterpolate(cfg,sourcestm,mri_resliced);
sourcetoneInt = ft_sourceinterpolate(cfg,sourcetone,mri_resliced);
cfg = [];
cfg.parameter= 'pow';
sourceDiffInt = ft_sourceinterpolate(cfg,sourceDiff,mri_resliced);
%% Visualize
insourceInt = sourceDiffInt;
maxval = max(insourceInt.pow)
minval = min(insourceInt.pow)
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.anaparameter = 'anatomy';
cfg.funcolorlim = [minval 0];
cfg.opacitylim = [minval 0];
cfg.opacitymap = 'rampup';
ft_sourceplot(cfg, insourceInt)

%% seedbased voxelwise
lh_seed_vox = [67   95   165]; %% left hemisphere
lh_seed_pos = ft_warp_apply(mri_resliced.transform, lh_seed_vox, 'homogeneous');
cfg = [];
cfg.funparameter = 'pow';
cfg.anaparameter = 'anatomy';
cfg.locationcoordinates = 'head';
cfg.location = lh_seed_pos;
ft_sourceplot(cfg, sourceDiffInt);

%% AAL
disp('realign AAL atlas')
%mni to ortho using ANTS
aalfname = '/group/prishah/work/tACS/refbrains/ROI_MNI_V4_nocerel.nii';
aal = ft_read_atlas(aalfname);
system(['cp /group/prishah/work/tACS/refbrains/ROI_MNI_V4_nocerel.txt ' ftanatdir '/.']);
cd(ftanatdir)

if ~exist(fullfile(ftanatdir, 'brain2_aal_inv.nii'))
    system(['mv ' ftanatdir '/ROI_MNI_V4_nocerel.txt ' ftanatdir '/brain2_aal_inv.txt']);

    system(['WarpImageMultiTransform 3 ' aalfname ' brain2_aal_inv.nii -R segmented.nii --use-NN mni_to_ortho_SYNWarp.nii.gz mni_to_ortho_SYNAffine.txt'])
    system(['3dresample -master /group/prishah/work/tACS/refbrains/ROI_MNI_V4_nocerel.nii -input brain_aal_inv.nii -prefix temp_brain_aal_inv.nii'])
    system(['3dcalc -a temp_brain_aal_inv.nii -expr ' char(39) 'a*ispositive(9000-a)' char(39) ' -prefix brain2_aal_inv.nii'])
    system('rm temp_brain_aal_inv.nii')
end
naal = ft_read_atlas('brain2_aal_inv.nii');
naal = ft_convert_units(naal, 'cm'); % Use SI Units
naal.coordsys = 'neuromag';
cfg            = [];
cfg.interpmethod = 'nearest';      
cfg.parameter  = 'tissue';
aalInt  = ft_sourceinterpolate(cfg, naal, sourcemodel); 
%% Now perform LCMV to obtain trialwise virtual signals at each grid location

cfg = [];
cfg.covariance = 'yes';                
cfg.covariancewindow = 'all';
cfg.keeptrials = 'yes';
cdatastm = ft_timelockanalysis(cfg, ftdatastm);
cdatatone= ft_timelockanalysis(cfg, ftdatatone);

%% not based on single trial
cfg                   = [];
cfg.method            = 'lcmv';            
cfg.sourcemodel       = sourcemodel;
cfg.headmodel         = headmodel;
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.lambda       = '5%';
cfg.lcmv.keepfilter   = 'yes';                
cfg.lcmv.keepmom      = 'yes';
cfg.keeptrial         = 'yes';
cfg.lcmv.fixedori     = 'yes';
cfg.lcmv.projectmom   = 'yes';
csourcestm            = ft_sourceanalysis(cfg, cdatastm);
csourcetone           = ft_sourceanalysis(cfg, cdatatone);
save(fullfile(ftoutdir, 'source_pstm.mat'),'csourcestm')
save(fullfile(ftoutdir, 'source_tone.mat'),'csourcetone')

%% extract virtual time series

chansel  = ft_channelselection('MEG', ftdatastm.label); % find MEG sensor names
chanindx = match_str(ftdatastm.label, chansel);         % find MEG sensor indices
tic;
virtstm = cell(1,length(csourcestm.inside));
for kk = 1:length(csourcestm.inside)   
    if isempty(csourcestm.avg.filter{kk})
        continue;
    else
        for tt = 1:length(ftdatastm.trial)
            virtstm{kk}(tt,:) = csourcestm.avg.filter{kk} * ftdatastm.trial{tt}(chanindx,:);
        end
    end
end
toc
save(fullfile(ftoutdir, 'source_alltrials_pstm.mat'),'virtstm', '-v7.3')

%load(fullfile(ftoutdir, 'source_tone.mat'))
%load(fullfile(ftoutdir, 'all_dataclean_tone_pilot2.mat'));
%ftdatatone=alldata; clear alldata
chansel  = ft_channelselection('MEG', ftdatatone.label); % find MEG sensor names
chanindx = match_str(ftdatatone.label, chansel);         % find MEG sensor indices
tic;
virttone = cell(1,length(csourcetone.inside));
for kk = 1:length(csourcetone.inside)   
    if isempty(csourcetone.avg.filter{kk})
        continue;
    else,
        for tt = 1:length(ftdatatone.trial)
            virttone{kk}(tt,:) = csourcetone.avg.filter{kk} * ftdatatone.trial{tt}(chanindx,:);
        end
    end
end
toc
save(fullfile(ftoutdir, 'source_alltrials_tone.mat'),'virttone','-v7.3')
%% parcellation
cfg = [];
cfg.parcellation = 'tissue';
cfg.method = 'pca';
outputstm=[]; outputstm = ft_virtualchannel(cfg,cdatastm, csourcestm,aalInt);
outputtone=[]; outputtone = ft_virtualchannel(cfg,cdatatone, csourcetone,aalInt);

save(fullfile(ftoutdir, 'aal_virtch_pstm.mat'),'outputstm')
save(fullfile(ftoutdir, 'aal_virtch_tone.mat'),'outputtone')
%save outputstm; outputtone;
%%
cfg = [];
cfg.method = 'mtmfft';
cfg.trials = 'all';
cfg.output = 'fourier';%'powandcsd'
cfg.tapsmofrq = 5;
cfg.foilim = [fq fq];                  
freqstm=[]; freqstm = ft_freqanalysis(cfg, outputstm);
freqtone=[]; freqtone = ft_freqanalysis(cfg, outputtone);

cfg = [];
cfg.method = 'coh';
cfg.complex = 'absimag';
icohstm=ft_connectivityanalysis(cfg, freqstm);
icohtone = ft_connectivityanalysis(cfg, freqtone);

cfg = [];
cfg.method = 'wpli';
wplistm=ft_connectivityanalysis(cfg, freqstm);
wplitone=ft_connectivityanalysis(cfg, freqtone);

save(fullfile(ftoutdir, ['imagcoh_pstm_' num2str(fq) 'Hz_aal.mat']),'icohstm')
save(fullfile(ftoutdir, ['imagcoh_tone_' num2str(fq) 'Hz_aal.mat']),'icohtone')

save(fullfile(ftoutdir, ['wpli_pstm_' num2str(fq) 'Hz_aal.mat']),'wplistm')
save(fullfile(ftoutdir, ['wpli_tone_' num2str(fq) 'Hz_aal.mat']),'wplitone')