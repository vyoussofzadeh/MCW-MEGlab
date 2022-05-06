% clear, 
close all,
clc
% cd(outputdir)

%%
PN = load('./PN/par_meg.mat');
PN2 = load('./PN/PN_pow');
PN3 = load('./PN/PN_ROI');

DFN = load('./DFN/par_meg.mat');
DFN2 = load('./DFN/DFN_pow');
DFN3 = load('./DFN/DFN_ROI');

coor = PN.coor;

%%
outputdir_idivcompr = './indivcomp_parc';
if exist(outputdir_idivcompr, 'file') == 0, mkdir(outputdir_idivcompr); end
cd(outputdir_idivcompr)

%%
[C,ia,ib] = intersect(DFN2.Sub_all,PN2.Sub_all, 'stable');
disp(C')
a = [ia,ib];

par_DFN_sel  = DFN.par_indv(ia,:);
par_PN_sel = PN.par_indv(ib,:);

c = DFN2.pow(ia,:);
pow_PN_sel  = PN2.pow(ib,:);

L = size(pow_DFN_sel,1);

%% 
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
msk = 'pow';
projthresh = 0.6;

%%
addpath(allpath.connpath);
addpath(allpath.spm_path);

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

%% DFN
close all
% idx = [12:2:20,80:2:88];
% idx2 = idx-1;
% idx3 = [idx,idx2];

tmp = abs(mean(pow_DFN_sel,1));

tsk = 'dfn'; pow_sel = pow_DFN_sel;
D = PN2.source_diff_dics;
for i=1:size(pow_sel,1)
    
    tmp = abs(pow_sel(i,:)); 
%     tmp(tmp > 0) = NaN;
    tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
    D.(msk) = (tmp);
    %     cfg = [];
    %     cfg.mask = 'pow';
    %     cfg.loc = 'min';
    %     cfg.template = template_mri;
    %     cfg.savefile = [];
    %     cfg.volnorm     = 2; % yes: 1
    %     source_dics = vy_source_plot(cfg, D);
    [~, par_meg, ~] = vy_parcellate(D, atlas,'pow');
    %     par_meg.anatomy = par_meg.pow;
    %     vy_parcellate_plot(par_meg, coor, 'net');
    par_meg.pow (par_meg.pow < 0) = NaN;
    par_meg.pow (par_meg.pow < 0.7*max(par_meg.pow)) = NaN;
    par_meg1 = par_meg;
    
    %%
    [idx,ROI_peak] = max(abs(par_meg1.pow(idx2)));
    ROI_DFN(i) = ROI_peak;
    
    %%
%     par_meg1.pow = -par_meg1.pow;
    %     par_meg1.pow=zeros(116,1);
    %     par_meg1.pow(idx3)=2.*(par_meg.pow(idx3));
    
    %     source_dics1 = vy_vol_thresh(source_dics,projthresh,'pow'); % abs
    %     source_dics1.pow = -(source_dics1.pow./max(source_dics1.pow(:)));
    
    %%
    set(gcf,'Name',C{i}) %select the name you want
    savenname = [tsk,'_',num2str(i),'_',C{i}];
    vy_savenifti(par_meg1,'pow',[savenname,'_parc.nii']);
    %     print(['./', tsk,'/',names{i}],'-depsc');
    %     print(savenname,'-dpng');
    
    %-Surf-vis
    opt.run = [];
    opt.tsk = tsk;
    opt.subj = C{i};
    opt.savedir = [];
    opt.savenii = 2;
    %     opt.plot = '-mosaic';
    opt.plot = '-row';
    vy_surfce_vis([],[savenname,'_parc.nii'], opt);
    
    disp([num2str(i),'/', num2str(size(pow_sel,1))])
    pause(1)
    close all
end

%% PN
close all
% idx = [12:2:20,80:2:88];
% idx2 = idx-1;
% idx3 = [idx,idx2];

tsk = 'pn'; pow_sel = pow_PN_sel;
D = PN2.source_diff_dics;
for i=1:size(pow_sel,1)
    
    tmp = abs(pow_sel(i,:)); 
%     tmp(tmp > 0) = NaN;
    tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
    D.(msk) = (tmp);
    %     cfg = [];
    %     cfg.mask = 'pow';
    %     cfg.loc = 'min';
    %     cfg.template = template_mri;
    %     cfg.savefile = [];
    %     cfg.volnorm     = 2; % yes: 1
    %     source_dics = vy_source_plot(cfg, D);
    [~, par_meg, ~] = vy_parcellate(D, atlas,'pow');
    %     par_meg.anatomy = par_meg.pow;
    %     vy_parcellate_plot(par_meg, coor, 'net');
    par_meg.pow (par_meg.pow < 0) = NaN;
    par_meg.pow (par_meg.pow < 0.7*max(par_meg.pow)) = NaN;
    par_meg1 = par_meg;
    
    %%    
    [idx,ROI_peak] = max(abs(par_meg1.pow(idx2)));
    ROI_PN(i) = ROI_peak;
    
    %%
%     par_meg1.pow = -par_meg1.pow;
    %     par_meg1.pow=zeros(116,1);
    %     par_meg1.pow(idx3)=2.*(par_meg.pow(idx3));
    
    %     source_dics1 = vy_vol_thresh(source_dics,projthresh,'pow'); % abs
    %     source_dics1.pow = -(source_dics1.pow./max(source_dics1.pow(:)));
    
    set(gcf,'Name',C{i}) %select the name you want
    savenname = [tsk,'_',num2str(i),'_',C{i}];
    vy_savenifti(par_meg1,'pow',[savenname,'_parc.nii']);
    %     print(['./', tsk,'/',names{i}],'-depsc');
    %     print(savenname,'-dpng');
    
%     %-Surf-vis
%     opt.run = [];
%     opt.tsk = tsk;
%     opt.subj = C{i};
%     opt.savedir = [];
%     opt.savenii = 2;
%     %     opt.plot = '-mosaic';
%     opt.plot = '-row';
%     vy_surfce_vis([],[savenname,'_parc.nii'], opt);
%     
%     disp([num2str(i),'/', num2str(size(pow_sel,1))])
%     pause(1)
    close all
end
