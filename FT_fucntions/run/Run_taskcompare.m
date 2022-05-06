% clear, 
close all,
clc

%%
PN = load('./PN/par_meg.mat');
PN2 = load('./PN/PN_pow');
PN3 = load('./PN/PN_ROI');

DFN = load('./DFN/par_meg.mat');
DFN2 = load('./DFN/DFN_pow');
DFN3 = load('./DFN/DFN_ROI');

coor = PN.coor;

%%
cd(outputdir)
outputdir_tskcompr = './taskcomp';
if exist(outputdir_tskcompr, 'file') == 0, mkdir(outputdir_tskcompr); end
cd(outputdir_tskcompr)

%%
[C,ia,ib] = intersect(DFN2.Sub_all,PN2.Sub_all, 'stable');
disp(C')
a = [ia,ib];

par_DFN_sel  = DFN.par_indv(ia,:);
par_PN_sel = PN.par_indv(ib,:);

pow_DFN_sel = DFN2.pow(ia,:);
pow_PN_sel  = PN2.pow(ib,:);

%% 
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
msk = 'pow';

%% G-average
%--DFN
tsk = 'dfn';
msk = 'pow';

D = [];
D = DFN2.source_diff_dics;
mpow = squeeze(mean(atanh(pow_DFN_sel),1)); % fisher-score transformation
% mpow = vy_normalize(pow_sel);
D.(msk) = mpow';

%--
cfg = [];
cfg.mask = msk;
cfg.template = template_mri;
cfg.savefile = 'DFN_dics_group_1';
cfg.volnorm     = 2; % yes: 1
D_DFN = vy_source_plot(cfg, D);
set(gcf,'name','DFN','numbertitle','off')
savepath = ('groupave_DFN.mat');
save(savepath, 'D','mpow', '-v7.3');
savenii = 'DFN_pow.nii';
vy_savenifti(D_DFN, msk, savenii);

clear savepath
savepath{1} = 'DFN_dics_group_2';
savepath{2} = 'DFN_dics_group_3';
cfg = [];
cfg.subj = ['g_',tsk];
cfg.mask = 'pow';
cfg.thre = 0.6;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, D_DFN);

%--PN
tsk = 'pn';
D = PN2.source_diff_dics;
mpow = squeeze(mean(atanh(pow_PN_sel),1)); % fisher-score transformation
% mpow = vy_normalize(pow_sel);
D.(msk) = mpow';
%--
cfg = [];
cfg.mask = msk;
cfg.template = template_mri;
cfg.savefile = 'PN_dics_group_1';
cfg.volnorm     = 2; % yes: 1
D_PN = vy_source_plot(cfg, D);
set(gcf,'name','PN','numbertitle','off');
savepath = ('groupave_PN.mat');
save(savepath, 'D','mpow', '-v7.3');
savenii = 'PN_pow.nii';
vy_savenifti(D_PN, msk, savenii);

clear savepath
savepath{1} = 'PN_dics_group_2';
savepath{2} = 'PN_dics_group_3';
cfg = [];
cfg.subj = ['g_',tsk];
cfg.mask = msk;
cfg.thre = 0.6;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, D_PN);

%%
addpath(allpath.connpath);
addpath(allpath.spm_path);

%%
group_source_dfn = ft_read_mri('DFN_pow.nii');
group_source_pn = ft_read_mri('PN_pow.nii');

%%
projthresh = 0.70;
s_vol_dfn = vy_vol_thresh(group_source_dfn, projthresh, 'anatomy'); % abs
projthresh = 0.70;
s_vol_pn = vy_vol_thresh(group_source_pn, projthresh, 'anatomy'); % abs

Opt = [];
Opt.savenii = 1; Opt.savefig = 0;
Opt.savename = 'dfn_groupave_thre';
% Opt.view = '-mosaic';
Opt.view = '-row';
vy_surfce_vis2(s_vol_dfn,[Opt.savename,'.nii'], Opt);
conn_mesh_display([Opt.savename,'.nii'], '');
set(gcf,'menubar','figure');
view([-110,20])
pause, %print the data
% clear savefig, savefig('groupave_thre_dfn_left.fig')
view([110,20])
pause, %print the data
% clear savefig, savefig('groupave_thre_dfn_right.fig')

Opt = [];
Opt.savenii = 1; Opt.savefig = 0;
Opt.savename = 'pn_groupave_thre';
% Opt.view = '-mosaic';
Opt.view = '-row';
vy_surfce_vis2(s_vol_pn,[Opt.savename,'.nii'], Opt);
conn_mesh_display([Opt.savename,'.nii'], '');
set(gcf,'menubar','figure');
view([-110,20])
pause, %print the data
% clear savefig, savefig('groupave_thre_pn_left.fig')
view([110,20])
pause, %print the data
% clear savefig, savefig('groupave_thre_pn_right.fig')

%% Parcellation (G-average), method 1
% close all
% 
% tsk = 'DFN';
% msk = 'anatomy'; 
% [~, D_par_DFN_nii, coor] = vy_parcellate(group_source_dfn, atlas, msk); 
% D_par_DFN_nii.([msk,'dimord']) = 'chan';
% 
% savenii = [tsk,'_group_par_nii.nii'];
% vy_savenifti(D_par_DFN_nii, msk, savenii);
% 
% [ROI_DFN, ROI_sel] = vy_ROI_report(D_par_DFN_nii,.8, coor, msk);
% disp(ROI_sel)
% 
% %-
% tsk = 'PN';
% [~, D_par_PN_nii, coor] = vy_parcellate(group_source_pn, atlas, msk); 
% D_par_PN_nii.([msk,'dimord']) = 'chan';
% 
% savenii = [tsk,'_group_par_nii.nii'];
% vy_savenifti(D_par_PN_nii, msk, savenii);
% 
% [ROI_PN, ROI_sel] = vy_ROI_report(D_par_PN_nii,.8, coor, msk);
% disp(ROI_sel)
% 
% group_source_par_dfn = ft_read_mri('DFN_group_par_nii.nii');
% group_source_par_pn = ft_read_mri('PN_group_par_nii.nii');
% 
% projthresh = 0.85;
% s_vol_dfn = vy_vol_thresh(group_source_par_dfn, projthresh, 'anatomy'); % abs
% s_vol_pn = vy_vol_thresh(group_source_par_pn, projthresh, 'anatomy'); % abs
% 
% Opt = []; Opt.savenii  = 1; Opt.savefig = 0; Opt.savename = 'DFN_group_par_thre_nii'; Opt.view = '-mosaic';
% vy_surfce_vis2(s_vol_dfn,[Opt.savename,'.nii'], Opt);
% conn_mesh_display([Opt.savename,'.nii'], '');
% % set(gcf,'menubar','figure');
% view([-110,20]), % manually print as, 'DFN_par_left.jpg'
% pause, %print the data
% % clear savefig, savefig('groupave_thre_dfn_left.fig')
% view([110,20]), % manually print as, 'DFN_par_left.jpg'
% pause, %print the data
% % clear savefig, savefig('groupave_thre_dfn_right.fig')
% 
% Opt = []; Opt.savenii  = 1; Opt.savefig = 0; Opt.savename = 'PN_group_par_thre_nii'; Opt.view = '-mosaic';
% vy_surfce_vis2(s_vol_pn,[Opt.savename,'.nii'], Opt);
% conn_mesh_display([Opt.savename,'.nii'], '');
% % set(gcf,'menubar','figure');
% view([-110,20])
% pause, %print the data
% % clear savefig, savefig('groupave_thre_dfn_left.fig')
% view([110,20])
% pause, %print the data
% clear savefig, savefig('groupave_thre_dfn_right.fig')

%% Parcellation (G-average), method 2
% DFN
tsk = 'DFN';
msk = 'pow';
D = [];
D = DFN2.source_diff_dics;
mpow = squeeze(mean(atanh(pow_DFN_sel),1)); % fisher-score transformation
% mpow = vy_normalize(pow_sel);
D.(msk) = mpow';
[~, D_par_DFN, coor] = vy_parcellate(D, atlas, msk); D_par_DFN.powdimord = 'chan';

%-
savepath = fullfile(outputdir,[tsk,'_group_par.mat']);
save(savepath, 'D_par_DFN', '-v7.3');
load (savepath);

%-
% clear savepath
% savepath{1} = fullfile(outputdir,[tsk,'dics_par_1']);
% savepath{2} = fullfile(outputdir,[tsk,'dics_par_2']);
% 
% cfg = [];
% cfg.subj = [tsk,'_group'];
% cfg.mask = msk;
% cfg.thre = 0.8;
% cfg.savepath = savepath;
% vy_mapvisualisation(cfg, D_par_DFN);

%-
savenii = [tsk,'_group_par.nii'];
vy_savenifti(D_par_DFN, msk, savenii);

tmp = abs(D_par_DFN.pow);
tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
D_par_DFN.pow = tmp;
D_par_DFN.pow = 5.8.*D_par_DFN.pow;
[ROI_DFN, ROI_sel] = vy_ROI_report(D_par_DFN,.6, coor, msk);
disp(ROI_sel)
savepath = [tsk,'_ROIs'];
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

textfile_rej = [tsk,'_ROI_par'];
writetable(ROI_sel,textfile_rej,'Delimiter',' ');

% PN
tsk = 'PN';
% msk = 'anatomy'; [~, D_par_PN, coor] = vy_parcellate(group_source_pn, atlas, msk); 
% D_par_PN.([msk,'dimord']) = 'chan';

D = PN2.source_diff_dics;
mpow = squeeze(mean(atanh(pow_PN_sel),1)); % fisher-score transformation
% mpow = vy_normalize(pow_sel);
D.(msk) = mpow';
[~, D_par_PN, coor] = vy_parcellate(D, atlas, msk);
D_par_PN.powdimord = 'chan';

savepath = fullfile(outputdir,[tsk,'_group_par.mat']);
save(savepath, 'D_par_PN', '-v7.3');
load (savepath);

%-
% clear savepath
% savepath{1} = fullfile(outputdir,[tsk,'dics_par_1']);
% savepath{2} = fullfile(outputdir,[tsk,'dics_par_2']);
% 
% cfg = [];
% cfg.subj = [tsk,'_group'];
% cfg.mask = msk;
% cfg.thre = 0.8;
% cfg.savepath = savepath;
% vy_mapvisualisation(cfg, D_par_PN);

%-
savenii = [tsk,'_group_par.nii'];
vy_savenifti(D_par_PN, msk, savenii);

tmp = abs(D_par_PN.pow);
tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
D_par_PN.pow = tmp;
D_par_PN.pow = 5.6.*D_par_PN.pow;

[ROI_PN, ROI_sel] = vy_ROI_report(D_par_PN,.6, coor, msk);
disp(ROI_sel)
savepath = [tsk,'_ROIs'];
hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

textfile_rej = [tsk,'_ROI_par'];
writetable(ROI_sel,textfile_rej,'Delimiter',' ');

group_source_par_dfn = ft_read_mri('DFN_group_par.nii');
group_source_par_pn = ft_read_mri('PN_group_par.nii');

projthresh = 0.85;
s_vol_dfn = vy_vol_thresh(group_source_par_dfn, projthresh, 'anatomy'); % abs
% s_vol_dfn.anatomy = (s_vol_dfn.anatomy)./(max(s_vol_dfn.anatomy(:)));
projthresh = 0.85;
s_vol_pn = vy_vol_thresh(group_source_par_pn, projthresh, 'anatomy'); % abs
% s_vol_pn.anatomy = (s_vol_pn.anatomy)./(max(s_vol_pn.anatomy(:)));

Opt = []; Opt.savenii  = 1; Opt.savefig = 0; Opt.savename = 'DFN_group_par_thre'; Opt.view = '-mosaic';
Opt.view = '-row';
vy_surfce_vis2(s_vol_dfn,[Opt.savename,'.nii'], Opt);
conn_mesh_display([Opt.savename,'.nii'], '');
% set(gcf,'menubar','figure');
view([-110,20]), % manually print as, 'DFN_par_left.jpg'
pause, %print the data
% clear savefig, savefig('groupave_thre_dfn_left.fig')
view([110,20]), % manually print as, 'DFN_par_left.jpg'
pause, %print the data
% clear savefig, savefig('groupave_thre_dfn_right.fig')

% s_vol_pn.anatomy = abs(s_vol_pn.anatomy);
Opt = []; Opt.savenii  = 1; Opt.savefig = 0; Opt.savename = 'PN_group_par_thre'; 
% Opt.view = '-mosaic'; 
Opt.view = '-row';
vy_surfce_vis2(s_vol_pn,[Opt.savename,'.nii'], Opt);
conn_mesh_display([Opt.savename,'.nii'], '');
% set(gcf,'menubar','figure');
view([-110,20])
pause, %print the data
% clear savefig, savefig('groupave_thre_dfn_left.fig')
view([110,20])
pause, %print the data
clear savefig, savefig('groupave_thre_dfn_right.fig')

%%
tsk = 'DFN'; msk = 'pow';
savepath = fullfile(outputdir,[tsk,'_group_par.mat']);
load(savepath);
[ROI_DFN, ROI_sel] = vy_ROI_report(D_par_DFN,.8, coor, msk);
% disp(ROI_sel)

tsk = 'PN';
savepath = fullfile(outputdir,[tsk,'_group_par.mat']);
load(savepath);
[ROI_PN, ROI_sel] = vy_ROI_report(D_par_PN,.8, coor, msk);
disp(ROI_sel)

%% ROI summary
%- parcellation, g-average
ROI_DFN1 = (ROI_DFN.(msk)); 
ROI_DFN1 = (ROI_DFN1)./(max(ROI_DFN1(:)));
ROI_PN1 = (ROI_PN.(msk)); 
ROI_PN1 = (ROI_PN1)./(max(ROI_PN1(:)));

% ROI_PN1 = ROI_PN; ROI_PN1.(msk) = (ROI_PN1.(msk))./(max(ROI_PN1.(msk)(:)));

% DataArray = [par_eeg_indv_CRM_sel(1,idx2)',par_eeg_indv_VGA_sel(1,idx2)'];
clear DataArray
DataArray = ([ROI_DFN1, ROI_PN1]);
DataArray1 = abs(DataArray);

Colors = [217,217,255;89,89,89]/256;
Colors = [32,100,166;227,159,157]/256;
save('ROIs.mat','ROI_DFN1', 'ROI_PN1','DataArray1');
load('ROIs.mat');

%% 
idx = [];
idx.central = [2,20]; 
idx.frontal = 4:2:18; 
idx.subcor = [22:2:48,78]; 
idx.Occ = 50:2:56; 
idx.pari= 58:2:70; 
idx.temp = 80:2:90;

% idx = [idx_central,idx_frontal,idx_subcor,idx_Occ,idx_pari,idx_temp];
idx1         = [idx.frontal,idx.temp,idx.pari];
idx1         = [idx.frontal];
idx1         = [idx.frontal, idx.temp];
idx1         = [idx.frontal, idx.temp, idx.pari];
% idx1 = 2:2:90;
idx2        = [idx1-1, idx1];

%% ROI mean
D = DFN2.source_diff_dics;
mpow = squeeze(mean(atanh(pow_DFN_sel),1)); % fisher-score transformation
% mpow = vy_normalize(pow_sel);
D.(msk) = mpow';
[~, D_par_DFN_sub, ~] = vy_parcellate(D, atlas, msk); D_par_DFN_sub.powdimord = 'chan';
[ROI_DFN_sub, ROI_sel_DFN] = vy_ROI_report(D_par_DFN_sub,.8, coor, msk);

D = PN2.source_diff_dics;
mpow = squeeze(mean(atanh(pow_PN_sel),1)); % fisher-score transformation
% mpow = pow_PN_sel(sub,:); % fisher-score transformation
D.(msk) = mpow';
[~, D_par_PN_sub, coor] = vy_parcellate(D, atlas, msk); D_par_PN_sub.powdimord = 'chan';
[ROI_PN_sub, ROI_sel_PN] = vy_ROI_report(D_par_PN_sub,.8, coor, msk);

ROI_DFN1 = abs(ROI_DFN_sub.(msk)); 
ROI_DFN1 = (ROI_DFN1)./(max(ROI_DFN1(:)));
ROI_PN1 = abs(ROI_PN_sub.(msk)); 
ROI_PN1 = (ROI_PN1)./(max(ROI_PN1(:)));

clear DataArray1
DataArray = ([ROI_DFN1, ROI_PN1]);
DataArray1 = abs(DataArray);

% idx1 = 2:2:90;
% idx2        = [idx1-1, idx1];
ytag = 'MEG source power';
xtag = 'ROI';
Run_roi_summary
savefig('ROI_summary_PN_vs_DFN.fig');

r = corr2(DataArray1(:,1),DataArray1(:,2));
% r = corr2(DataArray2(:,1),DataArray2(:,2));

%% Mapping
D = DFN2.source_diff_dics;
mpow = squeeze(mean(atanh(pow_DFN_sel),1)); % fisher-score transformation
% mpow = vy_normalize(pow_sel);
D.(msk) = mpow';
cfg = [];
cfg.mask = 'pow';
cfg.loc = 'min';
cfg.template = template_mri;
cfg.savefile = [];
cfg.volnorm     = 2; % yes: 1
source_dics_DFN = vy_source_plot(cfg, D);

mpow = squeeze(mean(atanh(pow_PN_sel),1)); % fisher-score transformation
% mpow = vy_normalize(pow_sel);
D.(msk) = mpow';
cfg = [];
cfg.mask = 'pow';
cfg.loc = 'min';
cfg.template = template_mri;
cfg.savefile = [];
cfg.volnorm     = 2; % yes: 1
source_dics_PN = vy_source_plot(cfg, D);

%% Mapping parcells.
D_par_DFN_sub.pow = ROI_DFN1;
D_par_DFN_sub.pow (D_par_DFN_sub.pow < 0.8*max(D_par_DFN_sub.pow)) = NaN;
vy_savenifti(D_par_DFN_sub,'pow','DFN_parc.nii');
%     print(['./', tsk,'/',names{i}],'-depsc');
%     print(savenname,'-dpng');

%-Surf-vis
opt.run = [];
opt.tsk = 'DFN';
opt.subj = 'all';
opt.savedir = [];
opt.savenii = 2;
%     opt.plot = '-mosaic';
opt.plot = '-row';
vy_surfce_vis([],'DFN_parc.nii', opt);

%%
D_par_PN_sub.pow = ROI_PN1;
D_par_PN_sub.pow (D_par_PN_sub.pow < 0.8*max(D_par_PN_sub.pow)) = NaN;
vy_savenifti(D_par_PN_sub,'pow','PN_parc.nii');
%     print(['./', tsk,'/',names{i}],'-depsc');
%     print(savenname,'-dpng');

%-Surf-vis
opt.run = [];
opt.tsk = 'PN';
opt.subj = 'all';
opt.savedir = [];
opt.savenii = 2;
%     opt.plot = '-mosaic';
opt.plot = '-row';
vy_surfce_vis([],'PN_parc.nii', opt);

%% ROI summary individuals
sub = 2;

D = DFN2.source_diff_dics;
mpow = pow_DFN_sel(sub,:); % fisher-score transformation
% mpow = vy_normalize(pow_sel);
D.(msk) = mpow';
[~, D_par_DFN_sub, ~] = vy_parcellate(D, atlas, msk); D_par_DFN_sub.powdimord = 'chan';
[ROI_DFN_sub, ROI_sel] = vy_ROI_report(D_par_DFN_sub,.8, coor, msk);

D = PN2.source_diff_dics;
mpow = pow_PN_sel(sub,:); % fisher-score transformation
D.(msk) = mpow';
[~, D_par_PN_sub, coor] = vy_parcellate(D, atlas, msk); D_par_PN_sub.powdimord = 'chan';
[ROI_PN_sub, ROI_sel] = vy_ROI_report(D_par_PN_sub,.8, coor, msk);

ROI_DFN1 = abs(ROI_DFN_sub.(msk)); 
ROI_DFN1 = (ROI_DFN1)./(max(ROI_DFN1(:)));
ROI_PN1 = abs(ROI_PN_sub.(msk)); 
ROI_PN1 = (ROI_PN1)./(max(ROI_PN1(:)));

clear DataArray1
DataArray = ([ROI_DFN1, ROI_PN1]);
DataArray1 = abs(DataArray);

% idx1 = 2:2:90;
% idx2        = [idx1-1, idx1];

ytag = 'MEG source power (t-value)';
xtag = 'ROI';
Run_roi_summary

%%
msk = 'pow';
idx3 = setdiff(1:116,idx2);

%-
tsk = 'DFN';
idx1 = find(DataArray1(:,1) < 0.8.*mu);
D_par_DFN1 = D_par_DFN;
D_par_DFN1.pow(idx3) = 0;
D_par_DFN1.pow(idx1) = 0;
ROI_DFN.label(DataArray1(:,1) > 0.8.*mu)
savenii = [tsk,'_group_par_thre.nii'];
vy_savenifti(D_par_DFN1, msk, savenii);
conn_mesh_display(savenii, '');

%-
tsk = 'PN';
idx1 = find(DataArray1(:,2) < 0.8.*mu);
D_par_PN1 = D_par_PN;
D_par_PN1.pow(idx3) = 0;
D_par_PN1.pow(idx1) = 0;
idx4 = find(abs(D_par_PN1.pow) > 0);
ROI_PN.label(DataArray1(:,2) > 0.8.*mu)
ROI_PN.label(idx4)
D_par_PN1.pow = abs(D_par_PN1.pow)./max(abs(D_par_PN1.pow(:)));
savenii = [tsk,'_group_par_thre.nii'];
vy_savenifti(D_par_PN1, msk, savenii);
conn_mesh_display(savenii, '');
view([-110,20]), % manually print as, 'DFN_par_left.jpg'
view([110,20])

%-mask
D_par_PN1 = D_par_PN;
D_par_PN1.pow = zeros(116,1);
% D_par_PN1.pow(idx2) = 1;
D_par_PN1.pow(81) = 1;
savenii = 'mask.nii';
vy_savenifti(D_par_PN1, msk, savenii);
conn_mesh_display(savenii, '');
ROI_PN.label(81)
view([-110,20]), % manually print as, 'DFN_par_left.jpg'

%%
%-PN
data = [];
data.value = abs(par_PN_sel);
data.label = PN.par_meg.label;
LI_PN = vy_laterality(data,idx2);

%-DFN
data = [];
data.value = abs(par_DFN_sel);
data.label = DFN.par_meg.label;
LI_DFN = vy_laterality(data, idx2);

%% BEST ROI for LATERALITY
k = 1; r_all =[];

for j=1:2:90
    
    idx2 = [j, j+1];
    
    %-PN
    data = [];
    data.value = abs(par_PN_sel);
    data.label = PN.par_meg.label;
    LI_PN = vy_laterality2(data,idx2);
    
    %-DFN
    data = [];
    data.value = abs(par_DFN_sel);
    data.label = DFN.par_meg.label;
    LI_DFN = vy_laterality2(data, idx2);
    
    Label{k} = DFN.par_meg.label{j};
    
    [r, p] = corr(LI_PN, LI_DFN);
    ml = mean((LI_DFN + LI_PN)/2)
    r_all(k) = r;
    ml_all(k) = ml;
    p_all(k) = p;
    k=1+k;
        
end

[val, idx] = sort(ml_all, 'descend');
disp(Label(idx(1:5)));
disp(r_all(idx(1:5)));
disp(p_all(idx(1:5)));
disp(ml_all(idx(1:5)));

idx3 = idx(1:5)*2-1;
idx2 = idx3+1;

%-PN
data = [];
data.value = abs(par_PN_sel);
data.label = PN.par_meg.label;
LI_PN = vy_laterality(data,idx2);

%-DFN
data = [];
data.value = abs(par_DFN_sel);
data.label = DFN.par_meg.label;
LI_DFN = vy_laterality(data, idx2);

r = corr(LI_PN, LI_DFN)

%% regression
mdl = fitlm(LI_PN, LI_DFN + 0.5.*LI_PN,'linear');
[r, p] = corr(LI_DFN + 0.5.*LI_PN, LI_DFN);

figure,
plot(mdl)
xlabel('Laterality');
set(gca,'color','none');
box off
set(gca,'FontName','HelveticaNeueLT Std Lt');
ylabel('Laterality','FontSize', 16);
xlabel('Subj','FontSize', 16);

ylabel([anal_tag, ' (mm)'])
title(['r=', num2str(corr(LI_DFN + 0.5.*LI_PN, LI_DFN)),' (p=', num2str(p),')']);


%% 
addpath('/data/MEG/Vahab/Github/MCW-MEGlab/tools/helpful tools/Violinplot-Matlab');
addpath('/data/MEG/Vahab/Github/MCW-MEGlab/tools/helpful tools/cbrewer');

% load carbig MPG Origin
% Origin = cellstr(Origin);
% figure
% vs = violinplot(MPG, Origin);
% % figure,
% % violinplot(DataArray)
% ylabel('Laterality','FontSize', 20);
% set(gca,'color','none');
% xlim([0.5, 2.5]);

%%
clear DataArray
DataArray = [LI_PN,LI_DFN];

Colors = [0.9 0 0;0 0.9 0];
% Colors = cbrewer('seq', 'PuBuGn', 2);

figure
UnivarScatter(DataArray,'Label',{'PN','DFN'},'MarkerFaceColor',Colors);
ylabel('Laterality','FontSize', 12);
set(gca,'color','none');

% figure,
% plotSpread(DataArray,'distributionMarkers',{'o'},'xNames', {'PN','DFN'},'categoryColors',{'r'});
% ylabel('Laterality','FontSize', 20);
% set(gcf, 'Position', [500   500   500  500]);
% set(gca,'color','none');

figure,bar(LI_PN);
view([90 -90])
L = length(LI_PN);
% set(gca,'Xtick', 1:L,'XtickLabel',C);
set(gca,'Xtick', 1:L,'XtickLabel',1:L);
% set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gcf, 'Position', [1500   500   500  500]);
grid on
set(gca,'color','none');
axis on
xlim([0,L+1])
xlabel('Subj');
ylabel('Laterality');
set(gcf, 'Position', [1500   500   500  500]);
title('PN')
clear savefig
savefig('LI_PN.fig')

figure,bar(LI_DFN);
view([90 -90])
L = length(LI_DFN);
% set(gca,'Xtick', 1:L,'XtickLabel',C);
set(gca,'Xtick', 1:L,'XtickLabel',1:L);
% set(gca,'FontSize',10,'XTickLabelRotation',90)
set(gcf, 'Position', [1500   500   500  500]);
grid on
set(gca,'color','none');
axis on
xlim([0,L+1])
xlabel('Subj');
ylabel('Laterality');
set(gcf, 'Position', [1500   500   500  500]);
title('DFN')
clear savefig
savefig('LI_DFN.fig')

%% Laterality index 
%-PN vs DFN
% close all; clc
figure,
% subplot 121
hbar = bar(DataArray);
view([90 -90])
L = length(DataArray);
for i=1:L
    S{i} = num2str(i);
end
set(gca,'Xtick',1:L,'XtickLabel',S);
xlim([0,L+1]);
ylim([-1,1]);
set(gca,'color','none');
box off
xlabel('Subj');
ylabel('Laterality');
% set(gcf, 'Position', [800   500   1000  500]);
% title('Word-recognition','fontsize',16)
legend([{'PN'};{'DFN'}]);
% set(gca,'FontName','Arial');
set(gca,'FontName','HelveticaNeueLT Std Lt');
ylabel('Laterality','FontSize', 16);
xlabel('Subj','FontSize', 16);
set(gcf, 'Position', [500   500  400  500]);

Colors = [0.4922    0.0039    0.9063;0.9922    0.8672    0.0039];
Colors = [0  0  1;0    1    0.9];
Colors = [217,217,255;89,89,89]/256;

figure
% subplot 122
[xPositions, yPositions, Label, RangeCut] = UnivarScatter(DataArray,'Label',{'PN','DFN'},'MarkerFaceColor',Colors);

ylabel('Laterality','FontSize', 16);
xlabel('Task','FontSize', 16);
set(gca,'color','none');
% title('Word-recognition','fontsize',16)
disp(['PN: ',num2str(mean(DataArray(:,1))), '+-', num2str(std(DataArray(:,1)))]);
disp(['DFN: ',num2str(mean(DataArray(:,2))), '+-', num2str(std(DataArray(:,2)))]);
set(gca,'FontName','HelveticaNeueLT Std Lt');
hold on
f = [xPositions, yPositions];
for j=1:length(f)
   line([f(j,1),f(j,2)],[f(j,3),f(j,4)]); 
end
ylim([-1 1])
for k = 1:numel(hbar)
set(hbar(k),'FaceColor',Colors(k,:))
end
set(gcf, 'Position', [500 500 400 500]);
% clear(savefig), savefig('')

for i=1:L
    Sub{i} = [num2str(i),'_',C{i}];
end
disp(Sub'),

%% Outliers
df = diff(DataArray')'; baddata = find(abs(df) > 0.5);
% baddata2 = find((DataArray(:,1).*DataArray(:,2)) < 0);
baddata = [3,7,13,32];
DataArray(baddata,:) = [];

%-
C1 = C;
for i=1:length(baddata)
    C1{baddata(i)} = '';
end
C1(strcmp('',C1)) = [];
save('baddata.mat','baddata','C1');

%%
clear Sub1
for i=1:length(C1)
    Sub1{i} = [num2str(i),'_',C1{i}];
end
disp(Sub1')

figure,
% subplot 121
hbar = bar(DataArray);
% view([90 -90])
L = length(DataArray);
clear S
for i=1:L
    S{i} = num2str(i);
end
set(gca,'Xtick',1:L,'XtickLabel',S);
xlim([0,L+1]);
ylim([-0.1,1]);
set(gca,'color','none');
box off
xlabel('Subj');
ylabel('Laterality');
% set(gcf, 'Position', [800   500   1000  500]);
% title('Word-recognition','fontsize',16)
legend([{'PN'};{'DFN'}]);
% set(gca,'FontName','Arial');
set(gca,'FontName','HelveticaNeueLT Std Lt');
ylabel('Laterality','FontSize', 16);
xlabel('Subject #','FontSize', 16);
set(gcf, 'Position', [500   500  400  500]);
Colors = [0.4922    0.0039    0.9063;0.9922    0.8672    0.0039];
Colors = [0  0  1;0    1    0.9];
Colors = [217,217,255;89,89,89]/256;

for k = 1:numel(hbar)
set(hbar(k),'FaceColor',Colors(k,:))
end
savefig('LI_group.fig')

figure
% subplot 122
[xPositions, yPositions, Label, RangeCut] = UnivarScatter(DataArray,'Label',{'PN','DFN'},'MarkerFaceColor',Colors);
ylabel('Laterality','FontSize', 16);
xlabel('Modality','FontSize', 16);
set(gca,'color','none');
% title('Word-recognition','fontsize',16)
disp(['PN: ',num2str(mean(DataArray(:,1))), '+-', num2str(std(DataArray(:,1)))]);
disp(['DFN: ',num2str(mean(DataArray(:,2))), '+-', num2str(std(DataArray(:,2)))]);
set(gca,'FontName','HelveticaNeueLT Std Lt');
hold on
f = [xPositions, yPositions];
for j=1:length(f)
   line([f(j,1),f(j,2)],[f(j,3),f(j,4)]); 
end
set(gcf, 'Position', [500 500 400 500]);
ylim([-0.2 1])
savefig('LI_group_2.fig');

%%
hold on
line([0 L], [0.25 0.25]);

%%
addpath(genpath('/data/MEG/Vahab/Github/MCW-MEGlab/tools/helpful tools/gramm'))
DataArray2 = [DataArray(:,1);DataArray(:,2)];
 
clear C1
X    = cell(1, size(DataArray2,1));
X(1:size(DataArray2,1)/2) = {'PN'};
X(size(DataArray2,1)/2+1:size(DataArray2,1)) = {'DFN'};

col = [2.*ones(1,L), 3.*ones(1,L)];

clear g
g(1,1) = gramm('x',X', 'y',DataArray2, 'color',col);

g(1,1).stat_boxplot();
% g(1,1).geom_jitter('width',0.4,'height',0);
g(1,1).set_title('laterality analysis');
g.set_names('column','Origin','x','Task','y','Laterality');

figure('Position',[500 500 400 500]);
g.draw();
savefig('LI_group_3.fig');

%% G-average
pow_DFN_sel1 = pow_DFN_sel;
pow_DFN_sel1(baddata,:) = [];

%--DFN
tsk = 'dfn';
msk = 'pow';
D = [];
D = DFN2.source_diff_dics;
mpow = squeeze(mean(atanh(pow_DFN_sel1),1)); % fisher-score transformation
% mpow = vy_normalize(pow_sel);
D.(msk) = mpow';
%--
cfg = [];
cfg.mask = msk;
cfg.template = template_mri;
cfg.savefile = 'DFN_dics_group_1';
cfg.volnorm     = 2; % yes: 1
D_DFN = vy_source_plot(cfg, D);
set(gcf,'name','DFN','numbertitle','off')
savepath = ('groupave_DFN.mat');
save(savepath, 'D','mpow', '-v7.3');
savenii = 'DFN_pow.nii';
vy_savenifti(D_DFN, msk, savenii);

clear savepath
savepath{1} = 'DFN_dics_group_2';
savepath{2} = 'DFN_dics_group_3';
cfg = [];
cfg.subj = ['g_',tsk];
cfg.mask = 'pow';
cfg.thre = 0.6;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, D_DFN);

%--PN
pow_PN_sel1 = pow_PN_sel;
pow_PN_sel1(baddata,:) = [];
tsk = 'pn';
D = PN2.source_diff_dics;
mpow = squeeze(mean(atanh(pow_PN_sel1),1)); % fisher-score transformation
% mpow = vy_normalize(pow_sel);
D.(msk) = mpow';
%--
cfg = [];
cfg.mask = msk;
cfg.template = template_mri;
cfg.savefile = 'PN_dics_group_1';
cfg.volnorm     = 2; % yes: 1
D_PN = vy_source_plot(cfg, D);
set(gcf,'name','PN','numbertitle','off');
savepath = ('groupave_PN.mat');
save(savepath, 'D','mpow', '-v7.3');
savenii = 'PN_pow.nii';
vy_savenifti(D_PN, msk, savenii);

clear savepath
savepath{1} = 'PN_dics_group_2';
savepath{2} = 'PN_dics_group_3';
cfg = [];
cfg.subj = ['g_',tsk];
cfg.mask = msk;
cfg.thre = 0.6;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, D_PN);

%% Visualise parcellation
% PN.par_indv(baddata,:) = [];
% addpath(allpath.spm_path);
% par_meg.anatomy = mean(PN.par_indv,1)';
% vy_parcellate_plot(par_meg, coor, 'PN_par');
% PN.roi = par_meg.anatomy(idx2);
% 
% DFN.par_indv(baddata,:) = [];
% par_meg.anatomy = mean(DFN.par_indv,1)';
% vy_parcellate_plot(par_meg, coor, 'DFN_par');
% DFN.roi = par_meg.anatomy(idx2);

%%
par_DFN_sel2 = par_DFN_sel;
par_PN_sel2  = par_PN_sel;

par_DFN_sel2(baddata,:)=[];
par_PN_sel2(baddata,:)=[];

DFN_sub = mean(par_DFN_sel2,1);
PN_sub = mean(par_PN_sel2,1);

% DFN_sub = mean(DFN.par_indv,1);
% PN_sub = mean(PN.par_indv,1);

DFN_sub = DFN_sub./max(DFN_sub);
PN_sub = PN_sub./max(PN_sub);

% DFN_sub = DFN_sub./max(DFN_sub);
% PN_sub = PN_sub./max(PN_sub);

clear DataArray
% DataArray = [par_eeg_indv_CRM_sel(1,idx2)',par_eeg_indv_VGA_sel(1,idx2)'];
DataArray = ([DFN_sub; PN_sub]');
DataArray1 = abs(DataArray);
Run_roi_summary
save('ROIs_sel.mat','DataArray2','label');
r = corr2(DataArray2(:,1),DataArray2(:,2));

%% Group
DFN_sub = mean(DFN.par_indv(:,idx2),1);
PN_sub = mean(PN.par_indv(:,idx2),1);

DFN_sub = DFN_sub./max(DFN_sub);
PN_sub = PN_sub./max(PN_sub);

DFN_sub = DFN_sub./max(DFN_sub);
PN_sub = PN_sub./max(PN_sub);

clear DataArray
% DataArray = [par_eeg_indv_CRM_sel(1,idx2)',par_eeg_indv_VGA_sel(1,idx2)'];
DataArray = ([DFN_sub; PN_sub]');
DataArray1 = abs(DataArray);

Colors = [217,217,255;89,89,89]/256;

figure(6)
bar_handle = bar(DataArray1,'grouped','BarWidth', 1);
set(bar_handle(1),'FaceColor',Colors(1,:))
set(bar_handle(2),'FaceColor',Colors(2,:))
L = length(DataArray1);
set(gca,'Xtick', 1:L,'XtickLabel',label); set(gca,'FontSize',10,'XTickLabelRotation',45)
% set(gca,'Xtick', 1:L,'XtickLabel',([1:10,1:10]));

r = corr2(DataArray1(:,1),DataArray1(:,2));

grid off
box off
set(gca,'color','none');
% axis on
xlim([0,L+1])
set(gca,'FontName','HelveticaNeueLT Std Lt');
xlabel('ROI','FontSize', 12);
ylabel('Source power','FontSize', 12);
set(gcf, 'Position', [1000   500   1000  400]);
% title('Group average')
legend({'PN','DFN'})
savefig('group_compare.fig')

coor1 = coor(idx2,:);

% [idx, l] = sort(DataArray1(:,1),'descend');
% label(l(1:5))
% 13*idx(1:5);
% coor1(l(1:5),:);
% 
% [idx, l] = sort(DataArray1(:,2),'descend');
% label(l(1:5))
% 13.*idx(1:5);
% coor1(l(1:5),:);

%%
addpath(allpath.connpath);
addpath(allpath.spm_path);
cd sub

%% individual
for Sub_sel=1:29
    close all
    %--DFN
    disp(C1{Sub_sel})
    
    %--
    tsk = 'DFN';
    D = DFN2.source_diff_dics;
    mpow = squeeze(mean(atanh(pow_DFN_sel1(Sub_sel,:)),1)); % fisher-score transformation
    D.(msk) = mpow';
    cfg = [];
    cfg.mask = 'pow';
    cfg.template = template_mri;
    cfg.savefile = [tsk,'_dics_group_1'];
    cfg.volnorm     = 2; % yes: 1
    D = vy_source_plot(cfg, D);
    set(gcf,'name',tsk,'numbertitle','off')
    savenii = ['s',num2str(Sub_sel),'_',tsk,'.nii'];
    vy_savenifti(D, 'pow', savenii);
    
    %--
    tsk = 'PN';
    D = PN2.source_diff_dics;
    mpow = squeeze(mean(atanh(pow_PN_sel1(Sub_sel,:)),1)); % fisher-score transformation
    D.(msk) = mpow';
    cfg = [];
    cfg.mask = 'pow';
    cfg.template = template_mri;
    cfg.savefile = [tsk,'_dics_group_1'];
    cfg.volnorm     = 2; % yes: 1
    D = vy_source_plot(cfg, D);
    set(gcf,'name',tsk,'numbertitle','off')
    savenii = ['s',num2str(Sub_sel),'_',tsk,'.nii'];
    vy_savenifti(D, 'pow', savenii);
    
    indv_source_dfn = ft_read_mri(['s',num2str(Sub_sel),'_DFN.nii']);
    indv_source_pn = ft_read_mri(['s',num2str(Sub_sel),'_PN.nii']);
    
    projthresh = 0.60;
    
    indv_source_dfn.anatomy = (indv_source_dfn.anatomy + 0.65.*group_source_dfn.anatomy )./2;
    s_vol_dfn = vy_vol_thresh(indv_source_dfn, projthresh, 'anatomy'); % abs
    
    indv_source_pn.anatomy = (indv_source_pn.anatomy + 0.65.*group_source_pn.anatomy)./2;
    s_vol_pn = vy_vol_thresh(indv_source_pn, projthresh, 'anatomy'); % abs
    
    Opt = [];
    Opt.savenii = 1; Opt.savefig = 0;
    Opt.view = '-mosaic'; Opt.view = '-row';
    Opt.savename = ['s',num2str(Sub_sel),'_DFN_thre'];
    % vy_surfce_vis2(s_vol_dfn,[Opt.savename,'.nii'], Opt);
    % conn_mesh_display([Opt.savename,'.nii'], '');
    vy_savenifti(s_vol_dfn, 'anatomy', [Opt.savename,'.nii']);
    
    Opt = [];
    Opt.savenii = 1; Opt.savefig = 0;
    Opt.view = '-mosaic'; Opt.view = '-row';
    Opt.savename = ['s',num2str(Sub_sel),'_PN_thre'];
    % vy_surfce_vis2(s_vol_pn,[Opt.savename,'.nii'], Opt);
    % conn_mesh_display([Opt.savename,'.nii'], '');
    vy_savenifti(s_vol_pn, 'anatomy', [Opt.savename,'.nii']);
    
    %-- overlaying maps
    disp(C1{Sub_sel})
    savenname_dfn = ['s',num2str(Sub_sel),'_DFN_thre.nii'];
    savenname_pn = ['s',num2str(Sub_sel),'_PN_thre.nii'];
    
    s_vol_dfn = ft_read_mri(savenname_dfn);
    s_vol_pn = ft_read_mri(savenname_pn);
    
    savenname_comp = ['comb_',num2str(Sub_sel),'_',C1{Sub_sel},'.nii'];
    spm_imcalc({savenname_dfn, savenname_pn},savenname_comp, 'i1.*(i1<0)-i2.*(i2<0)');
    
    %-Surf-vis
    opt.run = [];
    opt.tsk = 'comb';
    opt.subj = ['S',num2str(Sub_sel),'_',C{Sub_sel}];
    opt.savedir = [];
    opt.savenii = 2;
    opt.plot = '-row';
    vy_surfce_vis([],savenname_comp, opt);
    %     masked{k} = subj_fmri{j}; k=k+1;
end

%% Parcellation, individual
subj = 71;
DFN_sub = DFN.par_indv(subj,idx2)';
PN_sub = PN.par_indv(subj,idx2)';

clear DataArray
DataArray = ([DFN_sub,PN_sub]);
DataArray1 = abs(DataArray);

Colors = [217,217,255;89,89,89]/256;

figure
bar_handle = bar(DataArray1,'grouped','BarWidth', 1);
set(bar_handle(1),'FaceColor',Colors(1,:))
set(bar_handle(2),'FaceColor',Colors(2,:))
L = length(DataArray1);
set(gca,'Xtick', 1:L,'XtickLabel',label); set(gca,'FontSize',10,'XTickLabelRotation',45)
grid off
box off
set(gca,'color','none');
% axis on
xlim([0,L+1])
set(gca,'FontName','HelveticaNeueLT Std Lt');
xlabel('ROI','FontSize', 20);
ylabel('DICS-BF abs(power)','FontSize', 20);
set(gcf, 'Position', [1000   500   1000  400]);
title('Subject')
legend({'DFN','PN'})

%-DFN
data = [];
data.value = abs(DFN.par_indv(subj,:)); 
data.value = data.value./max(data.value);
data.label = DFN.par_meg.label;
LI_sub_DFN = vy_laterality(data);

%-PN
data = [];
data.value = abs(PN.par_indv(subj,:));
data.value = data.value./max(data.value);
data.label = PN.par_meg.label;
LI_sub_PN = vy_laterality(data);

DataArray = [LI_sub_DFN,LI_sub_PN];
figure,
subplot 221
bar(DataArray);
L = length(DataArray);
for i=1:L
    S{i} = num2str(i);
end
set(gca,'Xtick',1:L,'XtickLabel',S);
xlim([0,L+1]);
ylim([-1,1]);
set(gca,'color','none');
box off
xlabel('Subject #');
ylabel('Laterality');
% set(gcf, 'Position', [800   500   1000  500]);
% title('Word-recognition','fontsize',16)
legend([{'DFN'};{'PN'}]);
% set(gca,'FontName','Arial');
set(gca,'FontName','HelveticaNeueLT Std Lt');
ylabel('Laterality','FontSize', 16);
xlabel('Subj','FontSize', 16);
set(gcf, 'Position', [500   500  400  500]);

%%

