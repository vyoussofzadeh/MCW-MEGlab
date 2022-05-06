%%
outputdir = fullfile(outdir,'group_dics_2');
if exist(outputdir, 'file') == 0, mkdir(outputdir); end
cd(outputdir)

%%
PN = load('./PN/par_meg.mat');
PN2 = load('./PN/PN_pow');
PN3 = load('./PN/PN_ROI');

DFN = load('./DFN/par_meg.mat');
DFN2 = load('./DFN/DFN_pow');
DFN3 = load('./DFN/DFN_ROI');

coor = PN.coor;

%%
outputdir_stat = fullfile(outputdir,'stat');
if exist(outputdir_stat, 'file') == 0
    mkdir(outputdir_stat);   %create the directory
end
cd(outputdir_stat)

%%
% template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
load temp_grid % low-res
template_grid = ft_convert_units(template_grid, 'mm');

%% Combine subject-specific structures, group-stats
tsk_dfn = 'DFN';
% pow_sel = pow_DFN_sel;
pow_sel = DFN2.pow;
stat_group_DFN = vy_betweensub_stat(pow_sel,DFN2);
stat_group_DFN.pos     = template_grid.pos;
stat_group_DFN.dim     = template_grid.dim;
stat_group_DFN.inside  = template_grid.inside;

tsk_pn = 'PN';
pow_sel = PN2.pow;
stat_group_PN = vy_betweensub_stat(pow_sel,PN2);
stat_group_PN.pos     = template_grid.pos;
stat_group_PN.dim     = template_grid.dim;
stat_group_PN.inside  = template_grid.inside;

%% Group mean
pow = [];
D = source_dics;
mpow = squeeze(mean(atanh(pow_sel),1)); % fisher-score transformation
% mpow = vy_normalize(pow_sel);
D.(msk) = mpow';

%%
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %

savefig = [tsk_dfn,'_dics_stat_group_1'];
cfg = [];
cfg.mask = 'stat';
cfg.loc  = 'min';
cfg.template = template_mri;
cfg.savefile = savefig;
cfg.volnorm  = 2; % yes: 1
source_stat_dfn = vy_source_plot(cfg, stat_group_DFN);
clear savefig
savefig('group_stat_dfn.fig')

savefig = [tsk_pn,'_dics_stat_group_1'];
cfg = [];
cfg.mask = 'stat';
cfg.loc  = 'min';
cfg.template = template_mri;
cfg.savefile = savefig;
cfg.volnorm  = 2; % yes: 1
source_stat_pn = vy_source_plot(cfg, stat_group_PN);
clear savefig, savefig('group_stat_pn.fig')

%%
clear savepath
savepath{1} = [tsk_dfn,'_dics_stat_group_2'];
savepath{2} = [tsk_dfn,'_dics_stat_group_3'];

cfg = [];
cfg.subj = 'group';
cfg.mask = 'stat';
cfg.thre = 0.90;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, source_stat_dfn);

clear savepath
savepath{1} = [tsk_pn,'_dics_stat_group_2'];
savepath{2} = [tsk_pn,'_dics_stat_group_3'];

cfg = [];
cfg.subj = 'group';
cfg.mask = 'stat';
cfg.thre = 0.70;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, source_stat_pn);

%%
savenii_dnf = 'groupave_dfn.nii';
vy_savenifti(source_stat_dfn,'stat',savenii_dnf);

savenii_pn = 'groupave_pn.nii';
vy_savenifti(source_stat_pn,'stat',savenii_pn);

%%
addpath(allpath.connpath);
addpath(allpath.spm_path);

%%
group_source_dfn = ft_read_mri(savenii_dnf);
group_source_pn = ft_read_mri(savenii_pn);

% Opt = [];
% Opt.savenii = 0; Opt.savefig = 0;
% Opt.savename = [tsk,'_groupave1'];
% vy_surfce_vis2(group_source,[tsk,'_groupave.nii'], Opt);

%%
projthresh = 0.90;
s_vol_dfn = vy_vol_thresh(group_source_dfn, projthresh, 'anatomy'); % abs
projthresh = 0.90;
s_vol_pn = vy_vol_thresh(group_source_pn, projthresh, 'anatomy'); % abs

Opt = [];
Opt.savenii = 1; Opt.savefig = 0;
Opt.savename = [tsk_dfn,'_groupstat_thre'];
vy_surfce_vis2(s_vol_dfn,[tsk_dfn,'_groupstat_thre.nii'], Opt);
view([-110,20])
clear savefig, savefig('groupstat_thre_dfn_left.fig')
view([110,20])
clear savefig, savefig('groupstat_thre_dfn_right.fig')

Opt = [];
Opt.savenii = 1; Opt.savefig = 0;
Opt.savename = [tsk_pn,'_groupstat_thre'];
vy_surfce_vis2(s_vol_pn,[tsk_pn,'_groupstat_thre.nii'], Opt);
view([-110,20])
clear savefig, savefig('groupstat_thre_pn_left.fig')
view([110,20])
clear savefig, savefig('groupstat_thre_pn_right.fig')

%% mean,
group_source_mean = group_source_pn;
group_source_mean.anatomy = (group_source_pn.anatomy + group_source_dfn.anatomy)./2;

projthresh = 0.80;
s_vol_mean = vy_vol_thresh(group_source_mean, projthresh, 'anatomy'); % abs

Opt = [];
Opt.savenii = 1; Opt.savefig = 0;
Opt.savename = 'mean_groupstat_thre';
vy_surfce_vis2(s_vol_mean,[tsk_dfn,'_groupstat_thre.nii'], Opt);

